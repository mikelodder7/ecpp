//! Allocation-free ECPP support for fixed-width `crypto-bigint` integers.
//!
//! This module stores the certificate in a caller-selected, const-generic
//! number of slots. Neither allocation APIs nor heap-backed integer types are
//! used.

use core::array;

use crypto_bigint::{
    BitOps, CheckedAdd, CheckedSub, Limb, NonZero, Odd, RandomBits, RandomMod, Uint, WideWord,
    Word,
    modular::{FixedMontyForm, FixedMontyParams},
};
use rand_core::CryptoRng;

use crate::cm::{ClassPolynomial, DISCRIMINANTS};

/// Default maximum number of nodes in a fixed-capacity proof.
pub const DEFAULT_PROOF_NODES: usize = 64;

/// Errors from allocation-free proving and verification.
#[derive(Clone, Copy, Debug, Eq, PartialEq, thiserror::Error)]
pub enum Error {
    /// The input or requested capacity is invalid.
    #[error("invalid input: {0}")]
    InvalidInput(&'static str),
    /// A compositeness witness was found.
    #[error("candidate is composite")]
    Composite,
    /// The configured CM or factor search did not find another proof step.
    #[error("ECPP search exhausted")]
    SearchExhausted,
    /// The certificate does not have enough node slots.
    #[error("primality proof capacity exhausted")]
    CapacityExhausted,
    /// The supplied certificate is malformed or fails a required identity.
    #[error("invalid primality proof: {0}")]
    InvalidProof(&'static str),
}

/// Result type for allocation-free operations.
pub type Result<T> = core::result::Result<T, Error>;

/// Search limits used by the allocation-free prover.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProverOptions {
    /// The largest small prime stripped from candidate curve orders.
    pub trial_division_limit: u32,
    /// The number of random x-coordinates tried for each candidate curve.
    pub point_attempts: u32,
}

impl Default for ProverOptions {
    fn default() -> Self {
        Self {
            trial_division_limit: 10_000,
            point_attempts: 128,
        }
    }
}

/// A canonical integer held in a fixed-capacity big-endian byte array.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Natural<const BYTES: usize>([u8; BYTES]);

impl<const BYTES: usize> Natural<BYTES> {
    /// Creates a fixed-capacity value from a big-endian magnitude.
    pub fn from_be_slice(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > BYTES {
            return Err(Error::InvalidInput("integer exceeds certificate width"));
        }
        let mut output = [0u8; BYTES];
        output[BYTES - bytes.len()..].copy_from_slice(bytes);
        Ok(Self(output))
    }

    /// Returns the fixed-width big-endian encoding.
    pub const fn as_be_bytes(&self) -> &[u8; BYTES] {
        &self.0
    }
}

/// A short Weierstrass curve encoded independently of an arithmetic backend.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Curve<const BYTES: usize> {
    /// Linear coefficient.
    pub a: Natural<BYTES>,
    /// Constant coefficient.
    pub b: Natural<BYTES>,
}

/// An affine point encoded independently of an arithmetic backend.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Point<const BYTES: usize> {
    /// Affine x-coordinate.
    pub x: Natural<BYTES>,
    /// Affine y-coordinate.
    pub y: Natural<BYTES>,
}

/// One allocation-free elliptic Pocklington reduction from `n` to `q`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct EcppStep<const BYTES: usize> {
    /// Integer certified by this step, assuming `q` is prime.
    pub n: Natural<BYTES>,
    /// Curve over `Z/nZ`.
    pub curve: Curve<BYTES>,
    /// Point used by the elliptic Pocklington criterion.
    pub point: Point<BYTES>,
    /// Multiplier `m/q` used to produce a point of order divisible by `q`.
    pub cofactor: Natural<BYTES>,
    /// Prime certified by the following proof node.
    pub q: Natural<BYTES>,
}

/// A node in an allocation-free recursive certificate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProofNode<const BYTES: usize> {
    /// An elliptic-curve reduction.
    EllipticCurve(EcppStep<BYTES>),
    /// A deterministic 64-bit base case.
    SmallPrime(u64),
}

/// A stack-allocated primality proof with room for `NODES` nodes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PrimalityProof<const BYTES: usize, const NODES: usize = DEFAULT_PROOF_NODES> {
    nodes: [Option<ProofNode<BYTES>>; NODES],
    len: usize,
}

impl<const BYTES: usize, const NODES: usize> PrimalityProof<BYTES, NODES> {
    /// Returns the number of occupied certificate nodes.
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Reports whether the certificate has no nodes.
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns a certificate node by index.
    pub fn node(&self, index: usize) -> Option<&ProofNode<BYTES>> {
        self.nodes.get(index)?.as_ref()
    }

    /// Returns the integer certified by the first node.
    pub fn number(&self) -> Option<Natural<BYTES>> {
        match self.node(0)? {
            ProofNode::EllipticCurve(step) => Some(step.n),
            ProofNode::SmallPrime(prime) => Natural::from_be_slice(&prime.to_be_bytes()).ok(),
        }
    }

    /// Verifies the complete proof chain.
    pub fn verify_with<const LIMBS: usize>(&self) -> Result<()> {
        verify_proof::<LIMBS, BYTES, NODES>(self)
    }

    /// Verifies that this proof certifies `candidate`.
    pub fn verify_for<const LIMBS: usize>(&self, candidate: &Uint<LIMBS>) -> Result<()> {
        if self.number().as_ref() != Some(&encode_uint(candidate)?) {
            return Err(Error::InvalidProof("proof is for a different integer"));
        }
        self.verify_with::<LIMBS>()
    }

    /// Copies this stack proof into the backend-neutral heap proof format.
    ///
    /// The resulting certificate can be verified with any enabled heap-backed
    /// arithmetic backend, including boxed or fixed-width `crypto-bigint`,
    /// `num-bigint`, `rug`, and OpenSSL.
    #[cfg(feature = "alloc")]
    #[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
    pub fn to_alloc(&self) -> crate::PrimalityProof {
        let nodes = (0..self.len)
            .filter_map(|index| self.node(index))
            .map(|node| match node {
                ProofNode::SmallPrime(prime) => crate::ProofNode::SmallPrime(*prime),
                ProofNode::EllipticCurve(step) => {
                    crate::ProofNode::EllipticCurve(crate::EcppStep {
                        n: crate::Natural::from_be_bytes(step.n.as_be_bytes()),
                        curve: crate::Curve {
                            a: crate::Natural::from_be_bytes(step.curve.a.as_be_bytes()),
                            b: crate::Natural::from_be_bytes(step.curve.b.as_be_bytes()),
                        },
                        point: crate::Point {
                            x: crate::Natural::from_be_bytes(step.point.x.as_be_bytes()),
                            y: crate::Natural::from_be_bytes(step.point.y.as_be_bytes()),
                        },
                        cofactor: crate::Natural::from_be_bytes(step.cofactor.as_be_bytes()),
                        q: crate::Natural::from_be_bytes(step.q.as_be_bytes()),
                    })
                }
            })
            .collect();
        crate::PrimalityProof { nodes }
    }

    fn new() -> Self {
        Self {
            nodes: array::from_fn(|_| None),
            len: 0,
        }
    }

    fn push(&mut self, node: ProofNode<BYTES>) -> Result<()> {
        let slot = self
            .nodes
            .get_mut(self.len)
            .ok_or(Error::CapacityExhausted)?;
        *slot = Some(node);
        self.len += 1;
        Ok(())
    }
}

/// Constructs an allocation-free ECPP proof using `rng`.
pub fn prove_with_rng<
    const LIMBS: usize,
    const BYTES: usize,
    const NODES: usize,
    R: CryptoRng + ?Sized,
>(
    candidate: &Uint<LIMBS>,
    rng: &mut R,
) -> Result<PrimalityProof<BYTES, NODES>> {
    prove_with_options(candidate, rng, ProverOptions::default())
}

/// Constructs an allocation-free ECPP proof with explicit search limits.
pub fn prove_with_options<
    const LIMBS: usize,
    const BYTES: usize,
    const NODES: usize,
    R: CryptoRng + ?Sized,
>(
    candidate: &Uint<LIMBS>,
    rng: &mut R,
    options: ProverOptions,
) -> Result<PrimalityProof<BYTES, NODES>> {
    if LIMBS == 0 || NODES == 0 || *candidate < Uint::from(2u8) {
        return Err(Error::InvalidInput(
            "candidate and proof capacity must be nonzero",
        ));
    }
    if !is_probable_prime(candidate, options.trial_division_limit.max(53)) {
        return Err(Error::Composite);
    }

    let mut proof = PrimalityProof::new();
    let mut current = *candidate;
    loop {
        if let Some(small) = to_u64(&current) {
            if !is_prime_u64(small) {
                return Err(Error::Composite);
            }
            proof.push(ProofNode::SmallPrime(small))?;
            return Ok(proof);
        }
        if proof.len + 1 >= NODES {
            return Err(Error::CapacityExhausted);
        }
        let step = find_step(&current, rng, options)?;
        current = step.q;
        proof.push(ProofNode::EllipticCurve(encode_step(&step)?))?;
    }
}

/// Constructs a proof using the operating-system RNG.
#[cfg(feature = "getrandom")]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn prove<const LIMBS: usize, const BYTES: usize, const NODES: usize>(
    candidate: &Uint<LIMBS>,
) -> Result<PrimalityProof<BYTES, NODES>> {
    prove_with_rng(candidate, &mut rand_core::UnwrapErr(getrandom::SysRng))
}

/// Generates a fixed-width prime and its allocation-free proof using `rng`.
pub fn from_rng<
    const LIMBS: usize,
    const BYTES: usize,
    const NODES: usize,
    R: CryptoRng + ?Sized,
>(
    bit_length: u32,
    rng: &mut R,
) -> Result<crate::ProvedPrime<Uint<LIMBS>, PrimalityProof<BYTES, NODES>>> {
    if bit_length < 2 || bit_length > Uint::<LIMBS>::BITS {
        return Err(Error::InvalidInput(
            "bit length is outside the integer width",
        ));
    }
    loop {
        let mut candidate = Uint::random_bits(rng, bit_length);
        candidate.set_bit_vartime(bit_length - 1, true);
        candidate.set_bit_vartime(0, true);
        if !is_probable_prime(&candidate, 10_000) {
            continue;
        }
        match prove_with_rng(&candidate, rng) {
            Ok(proof) => {
                return Ok(crate::ProvedPrime {
                    prime: candidate,
                    proof,
                });
            }
            Err(Error::Composite | Error::SearchExhausted) => {}
            Err(error) => return Err(error),
        }
    }
}

/// Generates a fixed-width proved prime using the operating-system RNG.
#[cfg(feature = "getrandom")]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn new<const LIMBS: usize, const BYTES: usize, const NODES: usize>(
    bit_length: u32,
) -> Result<crate::ProvedPrime<Uint<LIMBS>, PrimalityProof<BYTES, NODES>>> {
    from_rng(bit_length, &mut rand_core::UnwrapErr(getrandom::SysRng))
}

/// Verifies an allocation-free proof against `candidate`.
pub fn verify<const LIMBS: usize, const BYTES: usize, const NODES: usize>(
    candidate: &Uint<LIMBS>,
    proof: &PrimalityProof<BYTES, NODES>,
) -> Result<()> {
    proof.verify_for(candidate)
}

/// Proves primality and returns `false` if construction fails.
pub fn check_with<
    const LIMBS: usize,
    const BYTES: usize,
    const NODES: usize,
    R: CryptoRng + ?Sized,
>(
    candidate: &Uint<LIMBS>,
    rng: &mut R,
) -> bool {
    prove_with_rng::<LIMBS, BYTES, NODES, _>(candidate, rng).is_ok()
}

#[derive(Clone, Copy)]
struct AffineCurve<const LIMBS: usize> {
    a: Uint<LIMBS>,
    b: Uint<LIMBS>,
}

#[derive(Clone, Copy)]
struct WorkingStep<const LIMBS: usize> {
    n: Uint<LIMBS>,
    curve: AffineCurve<LIMBS>,
    point: PointValue<LIMBS>,
    cofactor: Uint<LIMBS>,
    q: Uint<LIMBS>,
}

#[derive(Clone, Copy)]
struct PointValue<const LIMBS: usize> {
    x: Uint<LIMBS>,
    y: Uint<LIMBS>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum AffinePoint<const LIMBS: usize> {
    Infinity,
    Finite { x: Uint<LIMBS>, y: Uint<LIMBS> },
}

fn find_step<const LIMBS: usize, R: CryptoRng + ?Sized>(
    candidate: &Uint<LIMBS>,
    rng: &mut R,
    options: ProverOptions,
) -> Result<WorkingStep<LIMBS>> {
    if let Some(order) = candidate.checked_add(&Uint::ONE).into_option() {
        if rem_u32(candidate, 4) == 3 {
            let curve = AffineCurve {
                a: Uint::ONE,
                b: Uint::ZERO,
            };
            if let Some(step) = try_order(candidate, &order, &curve, None, rng, options)? {
                return Ok(step);
            }
        }
        if rem_u32(candidate, 3) == 2 {
            let curve = AffineCurve {
                a: Uint::ZERO,
                b: Uint::ONE,
            };
            if let Some(step) = try_order(candidate, &order, &curve, None, rng, options)? {
                return Ok(step);
            }
        }
    }

    for discriminant in DISCRIMINANTS {
        let Some((trace, _)) = cornacchia(candidate, discriminant.value) else {
            continue;
        };
        let Some(j_values) = j_invariants(candidate, discriminant.polynomial) else {
            continue;
        };
        for j in j_values.into_iter().flatten() {
            let base = curve_from_j(candidate, &j)?;
            let twist = quadratic_twist(candidate, &base)?;
            let Some(n_plus_one) = candidate.checked_add(&Uint::ONE).into_option() else {
                continue;
            };
            let orders = [
                n_plus_one.checked_sub(&trace).into_option(),
                n_plus_one.checked_add(&trace).into_option(),
            ];
            for order in orders.into_iter().flatten() {
                if let Some(step) = try_order(candidate, &order, &base, Some(&twist), rng, options)?
                {
                    return Ok(step);
                }
            }
        }
    }
    Err(Error::SearchExhausted)
}

fn try_order<const LIMBS: usize, R: CryptoRng + ?Sized>(
    candidate: &Uint<LIMBS>,
    order: &Uint<LIMBS>,
    curve: &AffineCurve<LIMBS>,
    twist: Option<&AffineCurve<LIMBS>>,
    rng: &mut R,
    options: ProverOptions,
) -> Result<Option<WorkingStep<LIMBS>>> {
    let Some(q) = split_order(candidate, order, options.trial_division_limit, rng) else {
        return Ok(None);
    };
    let Some(q_nonzero) = nonzero(q) else {
        return Ok(None);
    };
    let (cofactor, remainder) = order.div_rem_vartime(&q_nonzero);
    if remainder != Uint::ZERO {
        return Ok(None);
    }
    for candidate_curve in core::iter::once(curve).chain(twist) {
        if let Some(point) = find_point_of_order(
            candidate,
            candidate_curve,
            &cofactor,
            &q,
            rng,
            options.point_attempts,
        )? {
            let AffinePoint::Finite { x, y } = point else {
                return Ok(None);
            };
            return Ok(Some(WorkingStep {
                n: *candidate,
                curve: AffineCurve {
                    a: candidate_curve.a,
                    b: candidate_curve.b,
                },
                point: PointValue { x, y },
                cofactor,
                q,
            }));
        }
    }
    Ok(None)
}

fn split_order<const LIMBS: usize, R: CryptoRng + ?Sized>(
    candidate: &Uint<LIMBS>,
    order: &Uint<LIMBS>,
    trial_limit: u32,
    rng: &mut R,
) -> Option<Uint<LIMBS>> {
    let mut remaining = *order;
    for prime in 2..=trial_limit {
        if !is_prime_u32(prime) {
            continue;
        }
        while rem_u32(&remaining, prime) == 0 {
            remaining = div_u32(&remaining, prime);
        }
    }
    let fourth_root = ceil_fourth_root(candidate)?;
    let bound = fourth_root
        .checked_add(&Uint::ONE)
        .into_option()?
        .checked_mul(&fourth_root.checked_add(&Uint::ONE).into_option()?)
        .into_option()?;
    find_large_prime_factor(&remaining, candidate, &bound, trial_limit, rng, 0)
}

fn find_large_prime_factor<const LIMBS: usize, R: CryptoRng + ?Sized>(
    value: &Uint<LIMBS>,
    candidate: &Uint<LIMBS>,
    bound: &Uint<LIMBS>,
    trial_limit: u32,
    rng: &mut R,
    depth: usize,
) -> Option<Uint<LIMBS>> {
    if value <= bound || depth > 32 {
        return None;
    }
    if is_probable_prime(value, trial_limit.max(53)) {
        return (*value < *candidate).then_some(*value);
    }
    let factor = pollard_rho(value, rng)?;
    let factor_nonzero = nonzero(factor)?;
    let (other, remainder) = value.div_rem_vartime(&factor_nonzero);
    if remainder != Uint::ZERO {
        return None;
    }
    let (first, second) = if factor >= other {
        (factor, other)
    } else {
        (other, factor)
    };
    find_large_prime_factor(&first, candidate, bound, trial_limit, rng, depth + 1)
        .or_else(|| find_large_prime_factor(&second, candidate, bound, trial_limit, rng, depth + 1))
}

fn pollard_rho<const LIMBS: usize, R: CryptoRng + ?Sized>(
    value: &Uint<LIMBS>,
    rng: &mut R,
) -> Option<Uint<LIMBS>> {
    if !value.bit_vartime(0) {
        return Some(Uint::from(2u8));
    }
    let modulus = nonzero(*value)?;
    const MAX_ITERATIONS: usize = 250_000;
    for _ in 0..16 {
        let mut x = Uint::random_mod_vartime(rng, &modulus);
        let mut y = x;
        let constant = Uint::random_mod_vartime(rng, &modulus).add_mod(&Uint::ONE, &modulus);
        for _ in 0..MAX_ITERATIONS {
            x = rho_step(&x, &constant, &modulus);
            y = rho_step(&rho_step(&y, &constant, &modulus), &constant, &modulus);
            let difference = if x >= y {
                x.wrapping_sub(&y)
            } else {
                y.wrapping_sub(&x)
            };
            let divisor = difference.gcd_vartime(value);
            if divisor == Uint::ONE {
                continue;
            }
            if divisor != *value {
                return Some(divisor);
            }
            break;
        }
    }
    None
}

fn rho_step<const LIMBS: usize>(
    value: &Uint<LIMBS>,
    constant: &Uint<LIMBS>,
    modulus: &NonZero<Uint<LIMBS>>,
) -> Uint<LIMBS> {
    value
        .mul_mod_vartime(value, modulus)
        .add_mod(constant, modulus)
}

fn j_invariants<const LIMBS: usize>(
    candidate: &Uint<LIMBS>,
    polynomial: ClassPolynomial,
) -> Option<[Option<Uint<LIMBS>>; 2]> {
    match polynomial {
        ClassPolynomial::Linear(root) => Some([Some(mod_signed(root, candidate)?), None]),
        ClassPolynomial::Quadratic { constant, linear } => {
            let modulus = nonzero(*candidate)?;
            let linear_mod = mod_signed(linear, candidate)?;
            let constant_mod = mod_signed(constant, candidate)?;
            let four_constant = constant_mod.mul_mod_vartime(&Uint::from(4u8), &modulus);
            let discriminant = linear_mod
                .mul_mod_vartime(&linear_mod, &modulus)
                .sub_mod(&four_constant, &modulus);
            let square_root = modular_sqrt(&discriminant, candidate)?;
            let inverse_two = candidate
                .checked_add(&Uint::ONE)
                .into_option()?
                .shr_vartime(1);
            let minus_linear = if linear_mod == Uint::ZERO {
                Uint::ZERO
            } else {
                candidate.wrapping_sub(&linear_mod)
            };
            let first = minus_linear
                .sub_mod(&square_root, &modulus)
                .mul_mod_vartime(&inverse_two, &modulus);
            let second = minus_linear
                .add_mod(&square_root, &modulus)
                .mul_mod_vartime(&inverse_two, &modulus);
            Some([Some(first), Some(second)])
        }
    }
}

fn cornacchia<const LIMBS: usize>(
    candidate: &Uint<LIMBS>,
    discriminant: i16,
) -> Option<(Uint<LIMBS>, Uint<LIMBS>)> {
    let absolute = discriminant.unsigned_abs() as u32;
    let residue = mod_signed(discriminant as i128, candidate)?;
    let mut root = modular_sqrt(&residue, candidate)?;
    if root.bit_vartime(0) != (absolute & 1 == 1) {
        root = candidate.wrapping_sub(&root);
    }

    // Compute 2n mod root directly so the fixed-width temporary cannot overflow.
    let root_nonzero = nonzero(root)?;
    let n_mod_root = candidate.rem_vartime(&root_nonzero);
    let mut previous = root;
    let mut current = n_mod_root.add_mod(&n_mod_root, &root_nonzero);
    let sqrt_n = candidate.floor_sqrt_vartime();
    let square = sqrt_n.checked_mul(&sqrt_n).into_option()?;
    let remainder = candidate.checked_sub(&square).into_option()?;
    let mut limit = sqrt_n.checked_add(&sqrt_n).into_option()?;
    if remainder > sqrt_n {
        limit = limit.checked_add(&Uint::ONE).into_option()?;
    }
    while current > limit {
        let current_nonzero = nonzero(current)?;
        let next = previous.rem_vartime(&current_nonzero);
        previous = current;
        current = next;
    }

    let u_squared = current.checked_mul(&current).into_option()?;
    let v_squared = four_n_minus_square_div(candidate, &u_squared, absolute)?;
    let v = v_squared.checked_sqrt_vartime()?;
    Some((current, v))
}

fn four_n_minus_square_div<const LIMBS: usize>(
    n: &Uint<LIMBS>,
    square: &Uint<LIMBS>,
    divisor: u32,
) -> Option<Uint<LIMBS>> {
    let (low, high) = Uint::overflowing_shl_vartime_wide((*n, Uint::ZERO), 2)?;
    let (low, borrow) = low.borrowing_sub(square, Limb::ZERO);
    let (high, final_borrow) = high.borrowing_sub(&Uint::ZERO, borrow);
    if final_borrow != Limb::ZERO {
        return None;
    }
    let high_words = high.as_words();
    if high_words.iter().skip(1).any(|word| *word != 0) {
        return None;
    }
    let mut output = [0 as Word; LIMBS];
    let mut remainder = high_words.first().copied().unwrap_or(0) as WideWord;
    let divisor = divisor as WideWord;
    for index in (0..LIMBS).rev() {
        let wide = (remainder << Word::BITS) | low.as_words()[index] as WideWord;
        output[index] = (wide / divisor) as Word;
        remainder = wide % divisor;
    }
    (remainder == 0).then(|| Uint::from_words(output))
}

fn curve_from_j<const LIMBS: usize>(
    candidate: &Uint<LIMBS>,
    j: &Uint<LIMBS>,
) -> Result<AffineCurve<LIMBS>> {
    let modulus = nonzero(*candidate).ok_or(Error::Composite)?;
    let denominator = Uint::from(1728u16).sub_mod(j, &modulus);
    let inverse = denominator
        .invert_mod(&modulus)
        .into_option()
        .ok_or(Error::Composite)?;
    let k = j.mul_mod_vartime(&inverse, &modulus);
    Ok(AffineCurve {
        a: k.mul_mod_vartime(&Uint::from(3u8), &modulus),
        b: k.mul_mod_vartime(&Uint::from(2u8), &modulus),
    })
}

fn quadratic_twist<const LIMBS: usize>(
    candidate: &Uint<LIMBS>,
    curve: &AffineCurve<LIMBS>,
) -> Result<AffineCurve<LIMBS>> {
    let odd = Odd::new(*candidate).into_option().ok_or(Error::Composite)?;
    let modulus = odd.as_nz_ref();
    let mut non_residue = Uint::from(2u8);
    while jacobi(&non_residue, &odd) != -1 {
        non_residue = non_residue
            .checked_add(&Uint::ONE)
            .into_option()
            .ok_or(Error::Composite)?;
        if non_residue >= *candidate {
            return Err(Error::Composite);
        }
    }
    let square = non_residue.mul_mod_vartime(&non_residue, modulus);
    Ok(AffineCurve {
        a: curve.a.mul_mod_vartime(&square, modulus),
        b: curve
            .b
            .mul_mod_vartime(&square, modulus)
            .mul_mod_vartime(&non_residue, modulus),
    })
}

fn find_point_of_order<const LIMBS: usize, R: CryptoRng + ?Sized>(
    modulus: &Uint<LIMBS>,
    curve: &AffineCurve<LIMBS>,
    cofactor: &Uint<LIMBS>,
    q: &Uint<LIMBS>,
    rng: &mut R,
    attempts: u32,
) -> Result<Option<AffinePoint<LIMBS>>> {
    let modulus_nonzero = nonzero(*modulus).ok_or(Error::Composite)?;
    for _ in 0..attempts {
        let x = Uint::random_mod_vartime(rng, &modulus_nonzero);
        let rhs = x
            .mul_mod_vartime(&x, &modulus_nonzero)
            .mul_mod_vartime(&x, &modulus_nonzero)
            .add_mod(
                &curve.a.mul_mod_vartime(&x, &modulus_nonzero),
                &modulus_nonzero,
            )
            .add_mod(&curve.b, &modulus_nonzero);
        let Some(y) = modular_sqrt(&rhs, modulus) else {
            continue;
        };
        let point = AffinePoint::Finite { x, y };
        let q_point = scalar_mul(curve, modulus, cofactor, &point)?;
        if q_point == AffinePoint::Infinity {
            continue;
        }
        if scalar_mul(curve, modulus, q, &q_point)? == AffinePoint::Infinity {
            return Ok(Some(point));
        }
    }
    Ok(None)
}

fn point_add<const LIMBS: usize>(
    curve: &AffineCurve<LIMBS>,
    modulus: &Uint<LIMBS>,
    left: &AffinePoint<LIMBS>,
    right: &AffinePoint<LIMBS>,
) -> Result<AffinePoint<LIMBS>> {
    let (x1, y1) = match left {
        AffinePoint::Infinity => return Ok(*right),
        AffinePoint::Finite { x, y } => (x, y),
    };
    let (x2, y2) = match right {
        AffinePoint::Infinity => return Ok(*left),
        AffinePoint::Finite { x, y } => (x, y),
    };
    let modulus_nonzero = nonzero(*modulus).ok_or(Error::Composite)?;

    let slope = if x1 == x2 {
        if y1.add_mod(y2, &modulus_nonzero) == Uint::ZERO {
            return Ok(AffinePoint::Infinity);
        }
        if y1 != y2 {
            return Err(Error::Composite);
        }
        let denominator = y1.mul_mod_vartime(&Uint::from(2u8), &modulus_nonzero);
        let inverse = denominator
            .invert_mod(&modulus_nonzero)
            .into_option()
            .ok_or(Error::Composite)?;
        x1.mul_mod_vartime(x1, &modulus_nonzero)
            .mul_mod_vartime(&Uint::from(3u8), &modulus_nonzero)
            .add_mod(&curve.a, &modulus_nonzero)
            .mul_mod_vartime(&inverse, &modulus_nonzero)
    } else {
        let numerator = y2.sub_mod(y1, &modulus_nonzero);
        let denominator = x2.sub_mod(x1, &modulus_nonzero);
        let inverse = denominator
            .invert_mod(&modulus_nonzero)
            .into_option()
            .ok_or(Error::Composite)?;
        numerator.mul_mod_vartime(&inverse, &modulus_nonzero)
    };
    let x3 = slope
        .mul_mod_vartime(&slope, &modulus_nonzero)
        .sub_mod(x1, &modulus_nonzero)
        .sub_mod(x2, &modulus_nonzero);
    let y3 = slope
        .mul_mod_vartime(&x1.sub_mod(&x3, &modulus_nonzero), &modulus_nonzero)
        .sub_mod(y1, &modulus_nonzero);
    Ok(AffinePoint::Finite { x: x3, y: y3 })
}

fn scalar_mul<const LIMBS: usize>(
    curve: &AffineCurve<LIMBS>,
    modulus: &Uint<LIMBS>,
    scalar: &Uint<LIMBS>,
    point: &AffinePoint<LIMBS>,
) -> Result<AffinePoint<LIMBS>> {
    let mut output = AffinePoint::Infinity;
    let mut addend = *point;
    let mut scalar = *scalar;
    while scalar != Uint::ZERO {
        if scalar.bit_vartime(0) {
            output = point_add(curve, modulus, &output, &addend)?;
        }
        scalar = scalar.shr_vartime(1);
        if scalar != Uint::ZERO {
            addend = point_add(curve, modulus, &addend, &addend)?;
        }
    }
    Ok(output)
}

fn encode_step<const LIMBS: usize, const BYTES: usize>(
    step: &WorkingStep<LIMBS>,
) -> Result<EcppStep<BYTES>> {
    Ok(EcppStep {
        n: encode_uint(&step.n)?,
        curve: Curve {
            a: encode_uint(&step.curve.a)?,
            b: encode_uint(&step.curve.b)?,
        },
        point: Point {
            x: encode_uint(&step.point.x)?,
            y: encode_uint(&step.point.y)?,
        },
        cofactor: encode_uint(&step.cofactor)?,
        q: encode_uint(&step.q)?,
    })
}

fn decode_step<const LIMBS: usize, const BYTES: usize>(
    step: &EcppStep<BYTES>,
) -> Result<WorkingStep<LIMBS>> {
    Ok(WorkingStep {
        n: decode_uint(&step.n)?,
        curve: AffineCurve {
            a: decode_uint(&step.curve.a)?,
            b: decode_uint(&step.curve.b)?,
        },
        point: PointValue {
            x: decode_uint(&step.point.x)?,
            y: decode_uint(&step.point.y)?,
        },
        cofactor: decode_uint(&step.cofactor)?,
        q: decode_uint(&step.q)?,
    })
}

fn verify_proof<const LIMBS: usize, const BYTES: usize, const NODES: usize>(
    proof: &PrimalityProof<BYTES, NODES>,
) -> Result<()> {
    if proof.is_empty() {
        return Err(Error::InvalidProof("certificate is empty"));
    }
    let mut expected: Option<Uint<LIMBS>> = None;
    for index in 0..proof.len {
        let node = proof
            .node(index)
            .ok_or(Error::InvalidProof("certificate has an empty slot"))?;
        match node {
            ProofNode::SmallPrime(prime) => {
                if index + 1 != proof.len {
                    return Err(Error::InvalidProof("small-prime node must be last"));
                }
                if expected.is_some_and(|value| value != Uint::from(*prime)) {
                    return Err(Error::InvalidProof("certificate chain is disconnected"));
                }
                if !is_prime_u64(*prime) {
                    return Err(Error::InvalidProof("base case is not prime"));
                }
            }
            ProofNode::EllipticCurve(step) => {
                if index + 1 == proof.len {
                    return Err(Error::InvalidProof("certificate has no base case"));
                }
                let step = decode_step(step)?;
                verify_step(&step, expected.as_ref())?;
                expected = Some(step.q);
            }
        }
    }
    Ok(())
}

fn verify_step<const LIMBS: usize>(
    step: &WorkingStep<LIMBS>,
    expected: Option<&Uint<LIMBS>>,
) -> Result<()> {
    let n = step.n;
    if expected.is_some_and(|value| *value != n) {
        return Err(Error::InvalidProof("certificate chain is disconnected"));
    }
    if !n.bit_vartime(0) || n < Uint::from(3u8) || step.q >= n || step.q < Uint::from(2u8) {
        return Err(Error::InvalidProof("invalid ECPP step integers"));
    }
    let fourth_root = ceil_fourth_root(&n).ok_or(Error::InvalidProof("invalid proof bound"))?;
    let root_plus_one = fourth_root
        .checked_add(&Uint::ONE)
        .into_option()
        .ok_or(Error::InvalidProof("invalid proof bound"))?;
    let bound = root_plus_one
        .checked_mul(&root_plus_one)
        .into_option()
        .ok_or(Error::InvalidProof("invalid proof bound"))?;
    if step.q <= bound {
        return Err(Error::InvalidProof(
            "q is below the elliptic Pocklington bound",
        ));
    }
    if step.curve.a >= n || step.curve.b >= n || step.point.x >= n || step.point.y >= n {
        return Err(Error::InvalidProof("curve data is not reduced"));
    }
    let modulus = nonzero(n).ok_or(Error::InvalidProof("zero modulus"))?;
    let a_squared = step.curve.a.mul_mod_vartime(&step.curve.a, &modulus);
    let a_cubed = a_squared.mul_mod_vartime(&step.curve.a, &modulus);
    let b_squared = step.curve.b.mul_mod_vartime(&step.curve.b, &modulus);
    let discriminant = a_cubed.mul_mod_vartime(&Uint::from(4u8), &modulus).add_mod(
        &b_squared.mul_mod_vartime(&Uint::from(27u8), &modulus),
        &modulus,
    );
    if discriminant.gcd_vartime(&n) != Uint::ONE {
        return Err(Error::InvalidProof(
            "curve is singular modulo a divisor of n",
        ));
    }
    let rhs = step
        .point
        .x
        .mul_mod_vartime(&step.point.x, &modulus)
        .mul_mod_vartime(&step.point.x, &modulus)
        .add_mod(
            &step.curve.a.mul_mod_vartime(&step.point.x, &modulus),
            &modulus,
        )
        .add_mod(&step.curve.b, &modulus);
    if step.point.y.mul_mod_vartime(&step.point.y, &modulus) != rhs {
        return Err(Error::InvalidProof("certificate point is not on the curve"));
    }
    let curve = AffineCurve {
        a: step.curve.a,
        b: step.curve.b,
    };
    let point = AffinePoint::Finite {
        x: step.point.x,
        y: step.point.y,
    };
    let q_point = scalar_mul(&curve, &n, &step.cofactor, &point)
        .map_err(|_| Error::InvalidProof("point multiplication is not defined modulo n"))?;
    if q_point == AffinePoint::Infinity {
        return Err(Error::InvalidProof(
            "cofactor annihilates the certificate point",
        ));
    }
    if scalar_mul(&curve, &n, &step.q, &q_point)
        .map_err(|_| Error::InvalidProof("point multiplication is not defined modulo n"))?
        != AffinePoint::Infinity
    {
        return Err(Error::InvalidProof("q does not annihilate the point"));
    }
    Ok(())
}

/// Computes the Jacobi symbol `(value / modulus)` with the classic binary
/// algorithm.
///
/// `Uint::jacobi_symbol_vartime` in released crypto-bigint (through 0.7.5)
/// returns an incorrect sign for some inputs of four or more limbs
/// (RustCrypto/crypto-bigint#1295, fixed upstream but unreleased), so it
/// must not be used here.
fn jacobi<const LIMBS: usize>(value: &Uint<LIMBS>, modulus: &Odd<Uint<LIMBS>>) -> i8 {
    let mut value = value.rem_vartime(modulus.as_nz_ref());
    let mut modulus = *modulus.as_ref();
    let mut result = 1i8;
    while value != Uint::ZERO {
        while !value.bit_vartime(0) {
            value = value.shr_vartime(1);
            // For odd values, residue 3 or 5 modulo 8 is exactly when bits
            // one and two differ.
            if modulus.bit_vartime(1) != modulus.bit_vartime(2) {
                result = -result;
            }
        }
        core::mem::swap(&mut value, &mut modulus);
        // Both operands are odd here, so bit one alone decides whether each
        // is 3 modulo 4.
        if value.bit_vartime(1) && modulus.bit_vartime(1) {
            result = -result;
        }
        // The divisor is odd at every iteration, so it is never zero.
        let Some(modulus_nonzero) = nonzero(modulus) else {
            return 0;
        };
        value = value.rem_vartime(&modulus_nonzero);
    }
    if modulus == Uint::ONE { result } else { 0 }
}

fn modular_sqrt<const LIMBS: usize>(
    value: &Uint<LIMBS>,
    modulus: &Uint<LIMBS>,
) -> Option<Uint<LIMBS>> {
    if *value == Uint::ZERO {
        return Some(Uint::ZERO);
    }
    let odd = Odd::new(*modulus).into_option()?;
    if jacobi(value, &odd) != 1 {
        return None;
    }
    let params = FixedMontyParams::new_vartime(odd);
    let value_monty = FixedMontyForm::new(value, &params);
    if rem_u32(modulus, 4) == 3 {
        let exponent = modulus
            .checked_add(&Uint::ONE)
            .into_option()?
            .shr_vartime(2);
        return Some(value_monty.pow_vartime(&exponent).retrieve());
    }
    let mut odd_part = modulus.wrapping_sub(&Uint::ONE);
    let exponent = odd_part.trailing_zeros_vartime();
    odd_part = odd_part.shr_vartime(exponent);
    let mut non_residue = Uint::from(2u8);
    while jacobi(&non_residue, &odd) != -1 {
        non_residue = non_residue.checked_add(&Uint::ONE).into_option()?;
        if non_residue >= *modulus {
            return None;
        }
    }
    let mut c = FixedMontyForm::new(&non_residue, &params)
        .pow_vartime(&odd_part)
        .retrieve();
    let mut x = value_monty
        .pow_vartime(
            &odd_part
                .checked_add(&Uint::ONE)
                .into_option()?
                .shr_vartime(1),
        )
        .retrieve();
    let mut t = value_monty.pow_vartime(&odd_part).retrieve();
    let mut m = exponent;
    let modulus_nonzero = odd.as_nz_ref();
    while t != Uint::ONE {
        let mut i = 1u32;
        let mut power = t.mul_mod_vartime(&t, modulus_nonzero);
        while power != Uint::ONE {
            power = power.mul_mod_vartime(&power, modulus_nonzero);
            i += 1;
            if i >= m {
                return None;
            }
        }
        let exponent = Uint::<LIMBS>::ONE.shl_vartime(m - i - 1);
        let b = FixedMontyForm::new(&c, &params)
            .pow_vartime(&exponent)
            .retrieve();
        x = x.mul_mod_vartime(&b, modulus_nonzero);
        let b_squared = b.mul_mod_vartime(&b, modulus_nonzero);
        t = t.mul_mod_vartime(&b_squared, modulus_nonzero);
        c = b_squared;
        m = i;
    }
    Some(x)
}

fn is_probable_prime<const LIMBS: usize>(candidate: &Uint<LIMBS>, trial_limit: u32) -> bool {
    if *candidate < Uint::from(2u8) {
        return false;
    }
    for prime in 2..=trial_limit {
        if !is_prime_u32(prime) {
            continue;
        }
        if *candidate == Uint::from(prime) {
            return true;
        }
        if rem_u32(candidate, prime) == 0 {
            return false;
        }
    }
    const BASES: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    BASES.iter().all(|base| miller_rabin(candidate, *base))
}

fn miller_rabin<const LIMBS: usize>(candidate: &Uint<LIMBS>, base: u64) -> bool {
    let base = Uint::from(base);
    if base >= *candidate {
        return true;
    }
    let Some(odd_modulus) = Odd::new(*candidate).into_option() else {
        return false;
    };
    let params = FixedMontyParams::new_vartime(odd_modulus);
    let minus_one = candidate.wrapping_sub(&Uint::ONE);
    let exponent = minus_one.trailing_zeros_vartime();
    let odd = minus_one.shr_vartime(exponent);
    let mut value = FixedMontyForm::new(&base, &params)
        .pow_vartime(&odd)
        .retrieve();
    if value == Uint::ONE || value == minus_one {
        return true;
    }
    let modulus = odd_modulus.as_nz_ref();
    for _ in 1..exponent {
        value = value.mul_mod_vartime(&value, modulus);
        if value == minus_one {
            return true;
        }
        if value == Uint::ONE {
            return false;
        }
    }
    false
}

fn mod_signed<const LIMBS: usize>(value: i128, modulus: &Uint<LIMBS>) -> Option<Uint<LIMBS>> {
    if Uint::<LIMBS>::BITS < u128::BITS {
        return None;
    }
    let modulus_nonzero = nonzero(*modulus)?;
    let magnitude = Uint::<LIMBS>::from(value.unsigned_abs()).rem_vartime(&modulus_nonzero);
    if value >= 0 || magnitude == Uint::ZERO {
        Some(magnitude)
    } else {
        Some(modulus.wrapping_sub(&magnitude))
    }
}

fn ceil_fourth_root<const LIMBS: usize>(value: &Uint<LIMBS>) -> Option<Uint<LIMBS>> {
    let mut root = value.floor_sqrt_vartime().floor_sqrt_vartime();
    let square = root.checked_mul(&root).into_option()?;
    let fourth = square.checked_mul(&square).into_option()?;
    if fourth < *value {
        root = root.checked_add(&Uint::ONE).into_option()?;
    }
    Some(root)
}

fn encode_uint<const LIMBS: usize, const BYTES: usize>(
    value: &Uint<LIMBS>,
) -> Result<Natural<BYTES>> {
    if BYTES != Limb::BYTES * LIMBS {
        return Err(Error::InvalidInput(
            "certificate width does not match integer width",
        ));
    }
    Natural::from_be_slice(value.to_be_bytes().as_ref())
}

fn decode_uint<const LIMBS: usize, const BYTES: usize>(
    value: &Natural<BYTES>,
) -> Result<Uint<LIMBS>> {
    if BYTES != Limb::BYTES * LIMBS {
        return Err(Error::InvalidProof(
            "certificate width does not match integer width",
        ));
    }
    Ok(Uint::from_be_slice(value.as_be_bytes()))
}

fn nonzero<const LIMBS: usize>(value: Uint<LIMBS>) -> Option<NonZero<Uint<LIMBS>>> {
    NonZero::new(value).into_option()
}

fn rem_u32<const LIMBS: usize>(value: &Uint<LIMBS>, divisor: u32) -> u32 {
    let divisor = NonZero::new(Limb::from(divisor)).into_option();
    divisor.map_or(0, |divisor| value.rem_limb(divisor).0 as u32)
}

fn div_u32<const LIMBS: usize>(value: &Uint<LIMBS>, divisor: u32) -> Uint<LIMBS> {
    let Some(divisor) = NonZero::new(Limb::from(divisor)).into_option() else {
        return *value;
    };
    value.div_rem_limb(divisor).0
}

fn to_u64<const LIMBS: usize>(value: &Uint<LIMBS>) -> Option<u64> {
    if value.bits_vartime() > 64 {
        return None;
    }
    let mut output = 0u64;
    let bytes = value.to_be_bytes();
    for byte in bytes.as_ref().iter().rev().take(8).rev() {
        output = (output << 8) | u64::from(*byte);
    }
    Some(output)
}

fn is_prime_u32(candidate: u32) -> bool {
    if candidate < 2 {
        return false;
    }
    let mut divisor = 2u32;
    while divisor <= candidate / divisor {
        if candidate.is_multiple_of(divisor) {
            return false;
        }
        divisor += 1 + (divisor & 1);
    }
    true
}

fn is_prime_u64(candidate: u64) -> bool {
    if candidate < 2 {
        return false;
    }
    for prime in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if candidate == prime {
            return true;
        }
        if candidate.is_multiple_of(prime) {
            return false;
        }
    }
    [2u64, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022]
        .iter()
        .all(|base| miller_rabin_u64(candidate, *base % candidate))
}

fn miller_rabin_u64(candidate: u64, base: u64) -> bool {
    if base == 0 {
        return true;
    }
    let mut odd = candidate - 1;
    let exponent = odd.trailing_zeros();
    odd >>= exponent;
    let mut value = modpow_u64(base, odd, candidate);
    if value == 1 || value == candidate - 1 {
        return true;
    }
    for _ in 1..exponent {
        value = ((value as u128 * value as u128) % candidate as u128) as u64;
        if value == candidate - 1 {
            return true;
        }
    }
    false
}

fn modpow_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut output = 1u64;
    while exponent != 0 {
        if exponent & 1 == 1 {
            output = ((output as u128 * base as u128) % modulus as u128) as u64;
        }
        base = ((base as u128 * base as u128) % modulus as u128) as u64;
        exponent >>= 1;
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::U256;
    use rand::{SeedableRng, rngs::StdRng};

    #[test]
    fn allocation_free_base_case() {
        let mut rng = StdRng::seed_from_u64(41);
        let proof: PrimalityProof<{ U256::BYTES }, 8> =
            prove_with_rng(&U256::from(65_537u32), &mut rng).unwrap();
        proof.verify_for(&U256::from(65_537u32)).unwrap();

        let generated: crate::ProvedPrime<U256, PrimalityProof<{ U256::BYTES }, 8>> =
            from_rng(16, &mut rng).unwrap();
        generated.proof.verify_for(&generated.prime).unwrap();
    }

    #[test]
    fn jacobi_is_correct_where_jacobi_symbol_vartime_is_not() {
        use crypto_bigint::U512;

        // crypto-bigint 0.7.5's `jacobi_symbol_vartime` returns -1 for this
        // input; GMP, num-bigint, and the binary algorithm agree on +1.
        let mut value_bytes = [0u8; 64];
        value_bytes[54..]
            .copy_from_slice(&[0x60, 0x93, 0x39, 0x2a, 0xf9, 0x34, 0x22, 0x72, 0xbc, 0x3f]);
        let mut modulus_bytes = [0u8; 64];
        modulus_bytes[23..].copy_from_slice(&[
            0x28, 0x29, 0x56, 0x7a, 0x53, 0xf3, 0xea, 0x42, 0xc3, 0xe3, 0xd8, 0x35, 0x1a, 0x66,
            0xa7, 0x29, 0x0c, 0xd7, 0x52, 0x15, 0xd8, 0x57, 0x95, 0xd1, 0xc3, 0x02, 0x38, 0x20,
            0x9f, 0xbb, 0x15, 0x30, 0x42, 0xba, 0x4f, 0x2d, 0xf1, 0xc8, 0x8f, 0xf4, 0x47,
        ]);
        let value = U512::from_be_slice(&value_bytes);
        let modulus = Odd::new(U512::from_be_slice(&modulus_bytes)).unwrap();
        assert_eq!(jacobi(&value, &modulus), 1);

        let seven = Odd::new(U512::from(7u8)).unwrap();
        assert_eq!(jacobi(&U512::from(2u8), &seven), 1);
        assert_eq!(jacobi(&U512::from(3u8), &seven), -1);
        assert_eq!(jacobi(&U512::ZERO, &seven), 0);
        let one = Odd::new(U512::ONE).unwrap();
        assert_eq!(jacobi(&U512::from(5u8), &one), 1);
    }

    #[test]
    fn allocation_free_rejects_composite() {
        let mut rng = StdRng::seed_from_u64(43);
        let result: Result<PrimalityProof<{ U256::BYTES }, 8>> =
            prove_with_rng(&U256::from(561u32), &mut rng);
        assert_eq!(result, Err(Error::Composite));
    }

    #[test]
    fn allocation_free_multistep_proof() {
        let prime =
            U256::from_be_hex("00000000000000000000000000000000ffffffffffffffffffffffffffffff61");
        let mut rng = StdRng::seed_from_u64(47);
        let proof: PrimalityProof<{ U256::BYTES }, 64> = prove_with_rng(&prime, &mut rng).unwrap();
        assert!(proof.len() > 1);
        proof.verify_for(&prime).unwrap();

        #[cfg(feature = "alloc")]
        {
            let heap_proof = proof.to_alloc();
            heap_proof
                .verify_with::<crate::arithmetic::CryptoBigint>()
                .unwrap();
            #[cfg(feature = "num-bigint")]
            heap_proof
                .verify_with::<crate::arithmetic::NumBigint>()
                .unwrap();
            #[cfg(feature = "rug")]
            heap_proof.verify_with::<crate::arithmetic::Rug>().unwrap();
            #[cfg(feature = "openssl")]
            heap_proof
                .verify_with::<crate::arithmetic::OpenSsl>()
                .unwrap();
        }
    }
}
