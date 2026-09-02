//! Prime generation, proving, and certificate verification.

use alloc::vec;
use alloc::vec::Vec;

use num_bigint::BigUint;
use num_integer::Integer as _;
use num_traits::{One, ToPrimitive, Zero};
use rand_core::CryptoRng;

use crate::certificate::{Curve, EcppStep, Point, PrimalityProof, ProofNode};
use crate::cm::{DISCRIMINANTS, cornacchia, j_invariants};
use crate::math::{
    SMALL_PRIME_LIMIT, integer_sqrt, is_prime_u64, is_probable_prime, jacobi, mod_inverse, mod_sub,
    modular_sqrt, small_primes,
};
use crate::{Error, Integer, Natural, Result};

/// A generated prime together with its independently verifiable proof.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProvedPrime<T> {
    /// Prime value in the requested backend.
    pub prime: T,
    /// Atkin–Morain certificate for `prime`.
    pub proof: PrimalityProof,
}

impl<T> ProvedPrime<T> {
    /// Splits the generated value into its prime and proof.
    pub fn into_parts(self) -> (T, PrimalityProof) {
        (self.prime, self.proof)
    }
}

/// Bounds for certificate construction work.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProverOptions {
    /// Trial-division bound used to find the smooth cofactor of a curve order.
    pub trial_division_limit: u32,
    /// Random x-coordinates tried for each candidate curve.
    pub point_attempts: u32,
    /// Maximum number of ECPP reductions in one certificate.
    pub max_depth: usize,
}

impl Default for ProverOptions {
    fn default() -> Self {
        Self {
            trial_division_limit: SMALL_PRIME_LIMIT,
            point_attempts: 128,
            max_depth: 64,
        }
    }
}

/// Generates a proved prime using the operating system RNG.
#[cfg(feature = "getrandom")]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn new<T: Integer>(bit_length: usize) -> Result<ProvedPrime<T>> {
    from_rng(bit_length, &mut rand_core::UnwrapErr(getrandom::SysRng))
}

/// Generates a proved prime using `rng`.
pub fn from_rng<T: Integer, R: CryptoRng + ?Sized>(
    bit_length: usize,
    rng: &mut R,
) -> Result<ProvedPrime<T>> {
    if bit_length < 2 {
        return Err(Error::InvalidInput("bit length must be at least two"));
    }
    let byte_length = bit_length.div_ceil(8);
    let excess = byte_length * 8 - bit_length;
    let screening_primes = small_primes(SMALL_PRIME_LIMIT);
    loop {
        let mut bytes = vec![0u8; byte_length];
        rng.fill_bytes(&mut bytes);
        bytes[0] &= u8::MAX >> excess;
        bytes[0] |= 1u8 << (7 - excess);
        bytes[byte_length - 1] |= 1;
        let candidate = BigUint::from_bytes_be(&bytes);
        if !is_probable_prime(&candidate, &screening_primes) {
            continue;
        }
        let proof = match prove_biguint_with_options(&candidate, rng, ProverOptions::default()) {
            Ok(proof) => proof,
            Err(Error::SearchExhausted { .. } | Error::Composite) => continue,
            Err(error) => return Err(error),
        };
        let prime = T::from_be_bytes(&bytes).ok_or(Error::InvalidInput(
            "generated prime does not fit destination backend",
        ))?;
        return Ok(ProvedPrime { prime, proof });
    }
}

/// Constructs an Atkin–Morain proof using the operating system RNG.
#[cfg(feature = "getrandom")]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn prove<T: Integer>(candidate: &T) -> Result<PrimalityProof> {
    prove_with_rng(candidate, &mut rand_core::UnwrapErr(getrandom::SysRng))
}

/// Constructs an Atkin–Morain proof using `rng`.
pub fn prove_with_rng<T: Integer, R: CryptoRng + ?Sized>(
    candidate: &T,
    rng: &mut R,
) -> Result<PrimalityProof> {
    prove_with_options(candidate, rng, ProverOptions::default())
}

/// Constructs an Atkin–Morain proof with explicit search bounds.
pub fn prove_with_options<T: Integer, R: CryptoRng + ?Sized>(
    candidate: &T,
    rng: &mut R,
    options: ProverOptions,
) -> Result<PrimalityProof> {
    if candidate.is_negative() {
        return Err(Error::InvalidInput("candidate must be non-negative"));
    }
    let candidate = BigUint::from_bytes_be(&candidate.to_be_bytes());
    prove_biguint_with_options(&candidate, rng, options)
}

/// Constructs a proof from an unsigned big-endian magnitude using the OS RNG.
///
/// This is the dependency-free integration point for wrappers such as
/// `unknown_order::BigNumber`: pass the result of its `to_bytes()` method.
#[cfg(feature = "getrandom")]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn prove_be_bytes(candidate: &[u8]) -> Result<PrimalityProof> {
    prove(&Natural::from_be_bytes(candidate))
}

/// Verifies that `proof` certifies `candidate`.
pub fn verify<T: Integer>(candidate: &T, proof: &PrimalityProof) -> Result<()> {
    proof.verify_for(candidate)
}

/// Verifies a proof against an unsigned big-endian magnitude.
pub fn verify_be_bytes(candidate: &[u8], proof: &PrimalityProof) -> Result<()> {
    proof.verify_for(&Natural::from_be_bytes(candidate))
}

/// Proves primality and returns `false` on a composite or exhausted search.
#[cfg(feature = "getrandom")]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn check<T: Integer>(candidate: &T) -> bool {
    prove(candidate).is_ok()
}

/// RNG-explicit counterpart to [`check`].
pub fn check_with<T: Integer, R: CryptoRng + ?Sized>(candidate: &T, rng: &mut R) -> bool {
    prove_with_rng(candidate, rng).is_ok()
}

/// Alias for [`check`]; an ECPP proof is stronger than a probable-prime test.
#[cfg(feature = "getrandom")]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn strong_check<T: Integer>(candidate: &T) -> bool {
    check(candidate)
}

/// RNG-explicit counterpart to [`strong_check`].
pub fn strong_check_with<T: Integer, R: CryptoRng + ?Sized>(candidate: &T, rng: &mut R) -> bool {
    check_with(candidate, rng)
}

fn prove_biguint_with_options<R: CryptoRng + ?Sized>(
    candidate: &BigUint,
    rng: &mut R,
    options: ProverOptions,
) -> Result<PrimalityProof> {
    if candidate < &BigUint::from(2u8) {
        return Err(Error::Composite);
    }
    let primes = small_primes(options.trial_division_limit.max(53));
    if !is_probable_prime(candidate, &primes) {
        return Err(Error::Composite);
    }

    let mut nodes = Vec::new();
    let mut current = candidate.clone();
    for _ in 0..options.max_depth {
        if let Some(small) = current.to_u64() {
            if !is_prime_u64(small) {
                return Err(Error::Composite);
            }
            nodes.push(ProofNode::SmallPrime(small));
            return Ok(PrimalityProof { nodes });
        }
        let step = find_step(&current, rng, &primes, options.point_attempts)?;
        current = step.q.to_biguint();
        nodes.push(ProofNode::EllipticCurve(step));
    }
    Err(Error::SearchExhausted {
        candidate: Natural::from_biguint(&current),
    })
}

#[derive(Clone)]
struct AffineCurve {
    a: BigUint,
    b: BigUint,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum AffinePoint {
    Infinity,
    Finite { x: BigUint, y: BigUint },
}

fn find_step<R: CryptoRng + ?Sized>(
    candidate: &BigUint,
    rng: &mut R,
    primes: &[u32],
    point_attempts: u32,
) -> Result<EcppStep> {
    let mut search = StepSearch {
        candidate,
        rng,
        primes,
        point_attempts,
    };
    let order = candidate + 1u8;
    if candidate % 4u8 == BigUint::from(3u8) {
        let curve = AffineCurve {
            a: BigUint::one(),
            b: BigUint::zero(),
        };
        if let Some(step) = try_order(&mut search, &order, &curve, None)? {
            return Ok(step);
        }
    }
    if candidate % 3u8 == BigUint::from(2u8) {
        let curve = AffineCurve {
            a: BigUint::zero(),
            b: BigUint::one(),
        };
        if let Some(step) = try_order(&mut search, &order, &curve, None)? {
            return Ok(step);
        }
    }

    for discriminant in DISCRIMINANTS {
        let Some((trace, _)) = cornacchia(candidate, discriminant.value) else {
            continue;
        };
        let Some(j_invariants) = j_invariants(candidate, discriminant.polynomial) else {
            continue;
        };
        for j in j_invariants {
            let base = curve_from_j(candidate, &j)?;
            let twist = quadratic_twist(candidate, &base)?;
            for order in [candidate + 1u8 - &trace, candidate + 1u8 + &trace] {
                if let Some(step) = try_order(&mut search, &order, &base, Some(&twist))? {
                    return Ok(step);
                }
            }
        }
    }
    Err(Error::SearchExhausted {
        candidate: Natural::from_biguint(candidate),
    })
}

struct StepSearch<'a, R: ?Sized> {
    candidate: &'a BigUint,
    rng: &'a mut R,
    primes: &'a [u32],
    point_attempts: u32,
}

fn try_order<R: CryptoRng + ?Sized>(
    search: &mut StepSearch<'_, R>,
    order: &BigUint,
    curve: &AffineCurve,
    twist: Option<&AffineCurve>,
) -> Result<Option<EcppStep>> {
    let Some(q) = split_order(search.candidate, order, search.primes, search.rng) else {
        return Ok(None);
    };
    for candidate_curve in core::iter::once(curve).chain(twist) {
        if let Some(point) = find_point_of_order(
            search.candidate,
            candidate_curve,
            order,
            &q,
            search.rng,
            search.point_attempts,
        )? {
            return Ok(Some(EcppStep {
                n: Natural::from_biguint(search.candidate),
                curve: Curve {
                    a: Natural::from_biguint(&candidate_curve.a),
                    b: Natural::from_biguint(&candidate_curve.b),
                },
                point: match point {
                    AffinePoint::Finite { x, y } => Point {
                        x: Natural::from_biguint(&x),
                        y: Natural::from_biguint(&y),
                    },
                    AffinePoint::Infinity => return Ok(None),
                },
                order: Natural::from_biguint(order),
                q: Natural::from_biguint(&q),
            }));
        }
    }
    Ok(None)
}

fn split_order<R: CryptoRng + ?Sized>(
    candidate: &BigUint,
    order: &BigUint,
    primes: &[u32],
    rng: &mut R,
) -> Option<BigUint> {
    let mut remaining = order.clone();
    for &prime in primes {
        while &remaining % prime == BigUint::zero() {
            remaining /= prime;
        }
    }
    let fourth_root = ceil_fourth_root(candidate);
    let bound = (&fourth_root + 1u8).pow(2);
    find_large_prime_factor(&remaining, candidate, &bound, primes, rng, 0)
}

fn find_large_prime_factor<R: CryptoRng + ?Sized>(
    value: &BigUint,
    candidate: &BigUint,
    bound: &BigUint,
    primes: &[u32],
    rng: &mut R,
    depth: usize,
) -> Option<BigUint> {
    if value <= bound || depth > 32 {
        return None;
    }
    if is_probable_prime(value, primes) {
        return (value < candidate).then(|| value.clone());
    }
    let factor = pollard_rho(value, rng)?;
    let other = value / &factor;
    let (first, second) = if factor >= other {
        (factor, other)
    } else {
        (other, factor)
    };
    find_large_prime_factor(&first, candidate, bound, primes, rng, depth + 1)
        .or_else(|| find_large_prime_factor(&second, candidate, bound, primes, rng, depth + 1))
}

fn pollard_rho<R: CryptoRng + ?Sized>(value: &BigUint, rng: &mut R) -> Option<BigUint> {
    if value.is_even() {
        return Some(BigUint::from(2u8));
    }
    let one = BigUint::one();
    const BATCH: usize = 128;
    const MAX_ITERATIONS: usize = 1_000_000;
    for _ in 0..16 {
        let mut y = random_below(value, rng);
        let constant = random_below(value, rng) + 1u8;
        let mut divisor = one.clone();
        let mut power = 1usize;
        let mut iterations = 0usize;
        let mut x = y.clone();
        let mut saved_y = y.clone();
        while divisor == one && iterations < MAX_ITERATIONS {
            x.clone_from(&y);
            for _ in 0..power {
                y = ((&y * &y) + &constant) % value;
            }
            iterations += power;
            let mut offset = 0usize;
            while offset < power && divisor == one {
                saved_y.clone_from(&y);
                let count = BATCH.min(power - offset);
                let mut product = one.clone();
                for _ in 0..count {
                    y = ((&y * &y) + &constant) % value;
                    let difference = if x >= y { &x - &y } else { &y - &x };
                    product = (product * difference) % value;
                }
                divisor = product.gcd(value);
                offset += count;
                iterations += count;
            }
            power = power.saturating_mul(2);
        }
        if divisor == *value {
            divisor = one.clone();
            while divisor == one && iterations < MAX_ITERATIONS * 2 {
                saved_y = ((&saved_y * &saved_y) + &constant) % value;
                let difference = if x >= saved_y {
                    &x - &saved_y
                } else {
                    &saved_y - &x
                };
                divisor = difference.gcd(value);
                iterations += 1;
            }
        }
        if divisor != one && divisor != *value {
            return Some(divisor);
        }
    }
    None
}

fn ceil_fourth_root(value: &BigUint) -> BigUint {
    let mut root = integer_sqrt(&integer_sqrt(value));
    if root.pow(4) < *value {
        root += 1u8;
    }
    root
}

fn curve_from_j(candidate: &BigUint, j: &BigUint) -> Result<AffineCurve> {
    let denominator = mod_sub(&BigUint::from(1728u16), j, candidate);
    let inverse = mod_inverse(&denominator, candidate)?;
    let k = (j * inverse) % candidate;
    Ok(AffineCurve {
        a: (&k * 3u8) % candidate,
        b: (&k * 2u8) % candidate,
    })
}

fn quadratic_twist(candidate: &BigUint, curve: &AffineCurve) -> Result<AffineCurve> {
    let mut non_residue = BigUint::from(2u8);
    while jacobi(&non_residue, candidate) != -1 {
        non_residue += 1u8;
        if &non_residue >= candidate {
            return Err(Error::Composite);
        }
    }
    let square = (&non_residue * &non_residue) % candidate;
    Ok(AffineCurve {
        a: (&curve.a * &square) % candidate,
        b: (&curve.b * &square * non_residue) % candidate,
    })
}

fn find_point_of_order<R: CryptoRng + ?Sized>(
    modulus: &BigUint,
    curve: &AffineCurve,
    order: &BigUint,
    q: &BigUint,
    rng: &mut R,
    attempts: u32,
) -> Result<Option<AffinePoint>> {
    let cofactor = order / q;
    for _ in 0..attempts {
        let x = random_below(modulus, rng);
        let rhs = ((&x * &x % modulus) * &x + (&curve.a * &x) + &curve.b) % modulus;
        let Some(y) = modular_sqrt(&rhs, modulus) else {
            continue;
        };
        let point = AffinePoint::Finite { x, y };
        let q_point = scalar_mul(curve, modulus, &cofactor, &point)?;
        if q_point == AffinePoint::Infinity {
            continue;
        }
        if scalar_mul(curve, modulus, q, &q_point)? == AffinePoint::Infinity {
            return Ok(Some(point));
        }
    }
    Ok(None)
}

fn random_below<R: CryptoRng + ?Sized>(modulus: &BigUint, rng: &mut R) -> BigUint {
    let byte_length = modulus.bits().div_ceil(8) as usize;
    let excess = byte_length * 8 - modulus.bits() as usize;
    loop {
        let mut bytes = vec![0u8; byte_length];
        rng.fill_bytes(&mut bytes);
        bytes[0] &= u8::MAX >> excess;
        let value = BigUint::from_bytes_be(&bytes);
        if &value < modulus {
            return value;
        }
    }
}

fn point_add(
    curve: &AffineCurve,
    modulus: &BigUint,
    left: &AffinePoint,
    right: &AffinePoint,
) -> Result<AffinePoint> {
    let (x1, y1) = match left {
        AffinePoint::Infinity => return Ok(right.clone()),
        AffinePoint::Finite { x, y } => (x, y),
    };
    let (x2, y2) = match right {
        AffinePoint::Infinity => return Ok(left.clone()),
        AffinePoint::Finite { x, y } => (x, y),
    };

    let slope = if x1 == x2 {
        if (y1 + y2) % modulus == BigUint::zero() {
            return Ok(AffinePoint::Infinity);
        }
        if y1 != y2 {
            return Err(Error::Composite);
        }
        let denominator = (y1 * 2u8) % modulus;
        let inverse = mod_inverse(&denominator, modulus)?;
        (((x1 * x1 * 3u8) + &curve.a) * inverse) % modulus
    } else {
        let numerator = mod_sub(y2, y1, modulus);
        let denominator = mod_sub(x2, x1, modulus);
        (numerator * mod_inverse(&denominator, modulus)?) % modulus
    };
    let x3 = mod_sub(
        &mod_sub(&(&slope * &slope % modulus), x1, modulus),
        x2,
        modulus,
    );
    let y3 = mod_sub(&(&slope * mod_sub(x1, &x3, modulus) % modulus), y1, modulus);
    Ok(AffinePoint::Finite { x: x3, y: y3 })
}

fn scalar_mul(
    curve: &AffineCurve,
    modulus: &BigUint,
    scalar: &BigUint,
    point: &AffinePoint,
) -> Result<AffinePoint> {
    let mut output = AffinePoint::Infinity;
    let mut addend = point.clone();
    let mut scalar = scalar.clone();
    while !scalar.is_zero() {
        if scalar.is_odd() {
            output = point_add(curve, modulus, &output, &addend)?;
        }
        scalar >>= 1usize;
        if !scalar.is_zero() {
            addend = point_add(curve, modulus, &addend, &addend)?;
        }
    }
    Ok(output)
}

pub(crate) fn verify_proof(proof: &PrimalityProof) -> Result<()> {
    if proof.nodes.is_empty() {
        return Err(Error::InvalidProof("certificate is empty"));
    }
    let mut expected: Option<BigUint> = None;
    for (index, node) in proof.nodes.iter().enumerate() {
        match node {
            ProofNode::SmallPrime(prime) => {
                if index + 1 != proof.nodes.len() {
                    return Err(Error::InvalidProof("small-prime node must be last"));
                }
                if expected
                    .as_ref()
                    .is_some_and(|value| value != &BigUint::from(*prime))
                {
                    return Err(Error::InvalidProof("certificate chain is disconnected"));
                }
                if !is_prime_u64(*prime) {
                    return Err(Error::InvalidProof("base case is not prime"));
                }
            }
            ProofNode::EllipticCurve(step) => {
                if index + 1 == proof.nodes.len() {
                    return Err(Error::InvalidProof("certificate has no base case"));
                }
                verify_step(step, expected.as_ref())?;
                expected = Some(step.q.to_biguint());
            }
        }
    }
    Ok(())
}

fn verify_step(step: &EcppStep, expected: Option<&BigUint>) -> Result<()> {
    let n = step.n.to_biguint();
    let q = step.q.to_biguint();
    let order = step.order.to_biguint();
    if expected.is_some_and(|value| value != &n) {
        return Err(Error::InvalidProof("certificate chain is disconnected"));
    }
    if n.is_even() || n < BigUint::from(3u8) || q >= n || q < BigUint::from(2u8) {
        return Err(Error::InvalidProof("invalid ECPP step integers"));
    }
    if &order % &q != BigUint::zero() {
        return Err(Error::InvalidProof("q does not divide the curve order"));
    }
    let bound = (ceil_fourth_root(&n) + 1u8).pow(2);
    if q <= bound {
        return Err(Error::InvalidProof(
            "q is below the elliptic Pocklington bound",
        ));
    }

    let curve = AffineCurve {
        a: step.curve.a.to_biguint(),
        b: step.curve.b.to_biguint(),
    };
    let point = AffinePoint::Finite {
        x: step.point.x.to_biguint(),
        y: step.point.y.to_biguint(),
    };
    if curve.a >= n || curve.b >= n {
        return Err(Error::InvalidProof("curve coefficients are not reduced"));
    }
    let (x, y) = match &point {
        AffinePoint::Finite { x, y } if x < &n && y < &n => (x, y),
        _ => return Err(Error::InvalidProof("point coordinates are not reduced")),
    };
    let discriminant = ((curve.a.modpow(&BigUint::from(3u8), &n) * 4u8)
        + (curve.b.modpow(&BigUint::from(2u8), &n) * 27u8))
        % &n;
    if discriminant.gcd(&n) != BigUint::one() {
        return Err(Error::InvalidProof(
            "curve is singular modulo a divisor of n",
        ));
    }
    let rhs = ((x * x % &n) * x + (&curve.a * x) + &curve.b) % &n;
    if y * y % &n != rhs {
        return Err(Error::InvalidProof("certificate point is not on the curve"));
    }

    let cofactor = &order / &q;
    let q_point = scalar_mul(&curve, &n, &cofactor, &point)
        .map_err(|_| Error::InvalidProof("point multiplication is not defined modulo n"))?;
    if q_point == AffinePoint::Infinity {
        return Err(Error::InvalidProof(
            "cofactor annihilates the certificate point",
        ));
    }
    let result = scalar_mul(&curve, &n, &q, &q_point)
        .map_err(|_| Error::InvalidProof("point multiplication is not defined modulo n"))?;
    if result != AffinePoint::Infinity {
        return Err(Error::InvalidProof(
            "curve order does not annihilate the point",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{SeedableRng, rngs::StdRng};

    #[test]
    fn deterministic_u64_base_cases() {
        for prime in [2, 3, 5, 17, 65_537, 18_446_744_073_709_551_557] {
            assert!(is_prime_u64(prime));
        }
        for composite in [0, 1, 4, 9, 341, 561, u64::MAX] {
            assert!(!is_prime_u64(composite));
        }
    }

    #[test]
    fn prove_and_tamper_with_128_bit_prime() {
        let prime = BigUint::parse_bytes(b"340282366920938463463374607431768211297", 10).unwrap();
        let mut rng = StdRng::seed_from_u64(7);
        let proof = prove_with_rng(&prime, &mut rng).unwrap();
        proof.verify_for(&prime).unwrap();
        assert!(proof.nodes.len() >= 2);

        let mut changed = proof.clone();
        if let ProofNode::EllipticCurve(step) = &mut changed.nodes[0] {
            step.point.x = Natural::from_be_bytes(&[1]);
        }
        assert!(changed.verify().is_err());
    }

    #[test]
    #[ignore = "expensive fixed 256-bit ECPP regression"]
    fn prove_secp256k1_field_modulus() {
        let prime = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap();
        let mut rng = StdRng::seed_from_u64(23);
        let proof = prove_with_rng(&prime, &mut rng).unwrap();
        proof.verify_for(&prime).unwrap();
    }

    #[test]
    fn rejects_composites() {
        let mut rng = StdRng::seed_from_u64(11);
        for composite in [0u128, 1, 4, 9, 561, 3_215_031_751] {
            assert_eq!(prove_with_rng(&composite, &mut rng), Err(Error::Composite));
        }
    }

    #[test]
    fn rng_and_options_public_apis() {
        let mut rng = StdRng::seed_from_u64(13);
        let generated: ProvedPrime<u64> = from_rng(16, &mut rng).unwrap();
        let (prime, proof) = generated.into_parts();
        verify(&prime, &proof).unwrap();

        let options = ProverOptions {
            trial_division_limit: 100,
            point_attempts: 8,
            max_depth: 2,
        };
        let proof = prove_with_options(&65_537u64, &mut rng, options).unwrap();
        assert_eq!(
            proof.number(),
            Some(Natural::from_be_bytes(&65_537u64.to_be_bytes()))
        );
        assert!(check_with(&65_537u64, &mut rng));
        assert!(strong_check_with(&65_537u64, &mut rng));
        assert!(!check_with(&561u64, &mut rng));
    }

    #[cfg(feature = "getrandom")]
    #[test]
    fn operating_system_rng_public_apis() {
        let generated: ProvedPrime<u64> = new(16).unwrap();
        let (prime, proof) = generated.into_parts();
        proof.verify_for(&prime).unwrap();

        let proof = prove(&65_537u64).unwrap();
        verify(&65_537u64, &proof).unwrap();
        assert!(check(&65_537u64));
        assert!(strong_check(&65_537u64));
    }

    #[cfg(feature = "crypto-bigint")]
    #[test]
    fn crypto_bigint_uses_the_same_certificate() {
        use crypto_bigint::U256;

        let prime = U256::from_u128(65_537);
        let mut rng = StdRng::seed_from_u64(19);
        let proof = prove_with_rng(&prime, &mut rng).unwrap();
        proof.verify_for(&prime).unwrap();
    }

    #[cfg(feature = "rug")]
    #[test]
    fn rug_uses_the_same_certificate() {
        let prime = rug::Integer::from(65_537u32);
        let mut rng = StdRng::seed_from_u64(29);
        let proof = prove_with_rng(&prime, &mut rng).unwrap();
        proof.verify_for(&prime).unwrap();
    }

    #[cfg(feature = "openssl")]
    #[test]
    fn openssl_uses_the_same_certificate() {
        let prime = openssl::bn::BigNum::from_u32(65_537).unwrap();
        let mut rng = StdRng::seed_from_u64(30);
        let proof = prove_with_rng(&prime, &mut rng).unwrap();
        proof.verify_for(&prime).unwrap();
    }

    #[cfg(feature = "getrandom")]
    #[test]
    fn byte_api_supports_wrapper_backends() {
        let bytes = 65_537u32.to_be_bytes();
        let proof = prove_be_bytes(&bytes).unwrap();
        verify_be_bytes(&bytes, &proof).unwrap();
    }
}
