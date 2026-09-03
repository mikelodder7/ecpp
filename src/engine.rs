use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

use rand_core::CryptoRng;

use crate::arithmetic::{ArithmeticBackend, cmp_u64, from_u64, from_u128, to_u64};
use crate::cm::{ClassPolynomial, DISCRIMINANTS};
use crate::{Error, PrimalityProof, ProofNode, Result};

fn zero<B: ArithmeticBackend>() -> Result<B> {
    from_u64::<B>(0)
}

fn one<B: ArithmeticBackend>() -> Result<B> {
    from_u64::<B>(1)
}

fn add_u64<B: ArithmeticBackend>(value: &B, other: u64) -> Result<B> {
    value.clone() + &from_u64::<B>(other)?
}

fn mul_u64<B: ArithmeticBackend>(value: &B, other: u64) -> Result<B> {
    value.clone() * &from_u64::<B>(other)?
}

fn modular_sub<B: ArithmeticBackend>(left: &B, right: &B, modulus: &B) -> Result<B> {
    if left >= right {
        (left.clone() - right)? % modulus
    } else {
        let difference = ((right.clone() - left)? % modulus)?;
        if difference.is_zero() {
            zero::<B>()
        } else {
            modulus.clone() - &difference
        }
    }
}

fn modular_mul<B: ArithmeticBackend>(left: &B, right: &B, modulus: &B) -> Result<B> {
    (left.clone() * right)? % modulus
}

fn modular_signed<B: ArithmeticBackend>(value: i128, modulus: &B) -> Result<B> {
    let magnitude = (from_u128::<B>(value.unsigned_abs())? % modulus)?;
    if value >= 0 || magnitude.is_zero() {
        Ok(magnitude)
    } else {
        modulus.clone() - &magnitude
    }
}

fn is_square<B: ArithmeticBackend>(value: &B) -> Result<Option<B>> {
    let root = value.sqrt()?;
    Ok(((root.clone() * &root)? == *value).then_some(root))
}

fn modular_sqrt<B: ArithmeticBackend>(value: &B, modulus: &B) -> Result<Option<B>> {
    let value = (value.clone() % modulus)?;
    if value.is_zero() {
        return Ok(Some(zero::<B>()?));
    }
    if cmp_u64::<B>(modulus, 2)? == Ordering::Equal {
        return Ok(Some(value));
    }
    if value.jacobi(modulus)? != 1 {
        return Ok(None);
    }
    // A value is 3 modulo 4 exactly when its low two bits are set.
    if modulus.bit(0) && modulus.bit(1) {
        let exponent = (add_u64::<B>(modulus, 1)? >> 2)?;
        return Ok(Some(value.modular_pow(&exponent, modulus)?));
    }

    let one = one::<B>()?;
    let mut odd = (modulus.clone() - &one)?;
    let mut exponent = 0u32;
    while odd.is_even() {
        odd = (odd >> 1)?;
        exponent += 1;
    }
    let mut non_residue = from_u64::<B>(2)?;
    while non_residue.jacobi(modulus)? != -1 {
        non_residue = add_u64::<B>(&non_residue, 1)?;
        if &non_residue >= modulus {
            return Ok(None);
        }
    }
    let mut c = non_residue.modular_pow(&odd, modulus)?;
    let x_exponent = (add_u64::<B>(&odd, 1)? >> 1)?;
    let mut x = value.modular_pow(&x_exponent, modulus)?;
    let mut t = value.modular_pow(&odd, modulus)?;
    let mut m = exponent;
    while !t.is_one() {
        let mut i = 1u32;
        let mut power = modular_mul::<B>(&t, &t, modulus)?;
        while !power.is_one() {
            power = modular_mul::<B>(&power, &power, modulus)?;
            i += 1;
            if i >= m {
                return Ok(None);
            }
        }
        let two_power = (one.clone() << (m - i - 1) as usize)?;
        let b = c.modular_pow(&two_power, modulus)?;
        x = modular_mul::<B>(&x, &b, modulus)?;
        let b_squared = modular_mul::<B>(&b, &b, modulus)?;
        t = modular_mul::<B>(&t, &b_squared, modulus)?;
        c = b_squared;
        m = i;
    }
    Ok(Some(x))
}

fn small_primes(limit: u32) -> Vec<u32> {
    let mut composite = vec![false; limit as usize + 1];
    let mut primes = Vec::new();
    for candidate in 2..=limit {
        if !composite[candidate as usize] {
            primes.push(candidate);
            if candidate <= limit / candidate {
                let mut multiple = candidate * candidate;
                while multiple <= limit {
                    composite[multiple as usize] = true;
                    multiple += candidate;
                }
            }
        }
    }
    primes
}

fn is_probable_prime<B: ArithmeticBackend>(candidate: &B, primes: &[u32]) -> Result<bool> {
    if cmp_u64::<B>(candidate, 2)? == Ordering::Less {
        return Ok(false);
    }
    for &prime in primes {
        let prime_value = from_u64::<B>(u64::from(prime))?;
        if candidate == &prime_value {
            return Ok(true);
        }
        if (candidate.clone() % &prime_value)?.is_zero() {
            return Ok(false);
        }
    }
    if candidate.is_even() {
        return Ok(false);
    }
    for base in [
        2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
    ] {
        if !miller_rabin::<B>(candidate, &from_u64::<B>(base)?)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn miller_rabin<B: ArithmeticBackend>(candidate: &B, base: &B) -> Result<bool> {
    if base >= candidate {
        return Ok(true);
    }
    let one = one::<B>()?;
    let minus_one = (candidate.clone() - &one)?;
    let mut odd = minus_one.clone();
    let mut exponent = 0u32;
    while odd.is_even() {
        odd = (odd >> 1)?;
        exponent += 1;
    }
    let mut value = base.modular_pow(&odd, candidate)?;
    if value == one || value == minus_one {
        return Ok(true);
    }
    for _ in 1..exponent {
        value = modular_mul::<B>(&value, &value, candidate)?;
        if value == minus_one {
            return Ok(true);
        }
        if value == one {
            return Ok(false);
        }
    }
    Ok(false)
}

fn negate_mod<B: ArithmeticBackend>(value: &B, modulus: &B) -> Result<B> {
    if value.is_zero() {
        zero::<B>()
    } else {
        modulus.clone() - value
    }
}

fn quadratic_roots<B: ArithmeticBackend>(
    candidate: &B,
    linear_mod: &B,
    constant_mod: &B,
) -> Result<Option<(B, B)>> {
    let linear_squared = modular_mul::<B>(linear_mod, linear_mod, candidate)?;
    let four_constant = (mul_u64::<B>(constant_mod, 4)? % candidate)?;
    let discriminant = modular_sub::<B>(&linear_squared, &four_constant, candidate)?;
    let Some(square_root) = modular_sqrt::<B>(&discriminant, candidate)? else {
        return Ok(None);
    };
    let inverse_two = (add_u64::<B>(candidate, 1)? >> 1)?;
    let minus_linear = negate_mod::<B>(linear_mod, candidate)?;
    let first = modular_mul::<B>(
        &modular_sub::<B>(&minus_linear, &square_root, candidate)?,
        &inverse_two,
        candidate,
    )?;
    let second = modular_mul::<B>(&(minus_linear + &square_root)?, &inverse_two, candidate)?;
    Ok(Some((first, second)))
}

/// A monic cubic `x³ + quadratic·x² + linear·x + constant` with coefficients
/// reduced modulo the candidate, plus the precomputed reductions of `x³` and
/// `x⁴` used to multiply residue polynomials of degree below three.
struct CubicPolynomial<B: ArithmeticBackend> {
    quadratic: B,
    linear: B,
    constant: B,
    x_cubed: [B; 3],
    x_fourth: [B; 3],
}

impl<B: ArithmeticBackend> CubicPolynomial<B> {
    fn new(candidate: &B, quadratic: B, linear: B, constant: B) -> Result<Self> {
        let x_cubed = [
            negate_mod::<B>(&constant, candidate)?,
            negate_mod::<B>(&linear, candidate)?,
            negate_mod::<B>(&quadratic, candidate)?,
        ];
        // x⁴ ≡ (q² − l)·x² + (q·l − c)·x + q·c below the cubic.
        let q_squared = modular_mul::<B>(&quadratic, &quadratic, candidate)?;
        let x_fourth = [
            modular_mul::<B>(&quadratic, &constant, candidate)?,
            modular_sub::<B>(
                &modular_mul::<B>(&quadratic, &linear, candidate)?,
                &constant,
                candidate,
            )?,
            modular_sub::<B>(&q_squared, &linear, candidate)?,
        ];
        Ok(Self {
            quadratic,
            linear,
            constant,
            x_cubed,
            x_fourth,
        })
    }

    fn evaluate(&self, candidate: &B, x: &B) -> Result<B> {
        let mut accumulator = ((x.clone() + &self.quadratic)? % candidate)?;
        accumulator =
            ((modular_mul::<B>(&accumulator, x, candidate)? + &self.linear)? % candidate)?;
        (modular_mul::<B>(&accumulator, x, candidate)? + &self.constant)? % candidate
    }

    fn multiply(&self, candidate: &B, left: &[B; 3], right: &[B; 3]) -> Result<[B; 3]> {
        let mut raw = [
            zero::<B>()?,
            zero::<B>()?,
            zero::<B>()?,
            zero::<B>()?,
            zero::<B>()?,
        ];
        for (i, left_coefficient) in left.iter().enumerate() {
            for (j, right_coefficient) in right.iter().enumerate() {
                let product = modular_mul::<B>(left_coefficient, right_coefficient, candidate)?;
                raw[i + j] = ((raw[i + j].clone() + &product)? % candidate)?;
            }
        }
        let [d0, d1, d2, d3, d4] = raw;
        let mut output = [d0, d1, d2];
        for (index, coefficient) in output.iter_mut().enumerate() {
            let cubed = modular_mul::<B>(&d3, &self.x_cubed[index], candidate)?;
            let fourth = modular_mul::<B>(&d4, &self.x_fourth[index], candidate)?;
            *coefficient = (((coefficient.clone() + &cubed)? + &fourth)? % candidate)?;
        }
        Ok(output)
    }

    fn power(&self, candidate: &B, base: &[B; 3], exponent: &B) -> Result<[B; 3]> {
        let mut output = [one::<B>()?, zero::<B>()?, zero::<B>()?];
        let mut base = base.clone();
        let bits = exponent.bit_length();
        for index in 0..bits {
            if exponent.bit(index) {
                output = self.multiply(candidate, &output, &base)?;
            }
            if index + 1 < bits {
                base = self.multiply(candidate, &base, &base)?;
            }
        }
        Ok(output)
    }
}

fn trim_polynomial<B: ArithmeticBackend>(polynomial: &mut Vec<B>) {
    while polynomial
        .last()
        .is_some_and(|coefficient| coefficient.is_zero())
    {
        polynomial.pop();
    }
}

fn polynomial_rem<B: ArithmeticBackend>(
    candidate: &B,
    dividend: &[B],
    divisor: &[B],
) -> Result<Vec<B>> {
    let divisor_degree = divisor.len() - 1;
    let lead_inverse = divisor[divisor_degree].modular_inverse(candidate)?;
    let mut remainder = dividend.to_vec();
    while remainder.len() > divisor_degree {
        let top = remainder.len() - 1;
        let coefficient = remainder[top].clone();
        if !coefficient.is_zero() {
            let factor = modular_mul::<B>(&coefficient, &lead_inverse, candidate)?;
            let shift = top - divisor_degree;
            for (offset, divisor_coefficient) in divisor.iter().enumerate() {
                let subtrahend = modular_mul::<B>(&factor, divisor_coefficient, candidate)?;
                remainder[shift + offset] =
                    modular_sub::<B>(&remainder[shift + offset], &subtrahend, candidate)?;
            }
        }
        remainder.pop();
    }
    trim_polynomial::<B>(&mut remainder);
    Ok(remainder)
}

fn polynomial_gcd<B: ArithmeticBackend>(
    candidate: &B,
    mut left: Vec<B>,
    mut right: Vec<B>,
) -> Result<Vec<B>> {
    trim_polynomial::<B>(&mut left);
    trim_polynomial::<B>(&mut right);
    while !right.is_empty() {
        let remainder = polynomial_rem::<B>(candidate, &left, &right)?;
        left = right;
        right = remainder;
    }
    Ok(left)
}

/// Finds roots of a fully split monic cubic by Cantor–Zassenhaus splitting:
/// `gcd(H, (x + r)^((n−1)/2) − 1)` isolates the roots whose shifted values
/// are quadratic residues.
fn cubic_roots<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    candidate: &B,
    polynomial: &CubicPolynomial<B>,
    rng: &mut R,
) -> Result<Option<Vec<B>>> {
    let one = one::<B>()?;
    let exponent = ((candidate.clone() - &one)? >> 1)?;
    let cubic = vec![
        polynomial.constant.clone(),
        polynomial.linear.clone(),
        polynomial.quadratic.clone(),
        one.clone(),
    ];
    for _ in 0..16 {
        let shift = random_below::<B, R>(candidate, rng)?;
        let base = [shift, one.clone(), zero::<B>()?];
        let mut split = polynomial.power(candidate, &base, &exponent)?;
        split[0] = modular_sub::<B>(&split[0], &one, candidate)?;
        let factor = polynomial_gcd::<B>(candidate, cubic.clone(), split.to_vec())?;
        let candidates = match factor.len() {
            2 => {
                let root = modular_mul::<B>(
                    &negate_mod::<B>(&factor[0], candidate)?,
                    &factor[1].modular_inverse(candidate)?,
                    candidate,
                )?;
                // Deflate by the root: H = (x − r)(x² + p·x + q).
                let deflated_linear = ((polynomial.quadratic.clone() + &root)? % candidate)?;
                let deflated_constant = ((modular_mul::<B>(&root, &deflated_linear, candidate)?
                    + &polynomial.linear)?
                    % candidate)?;
                let mut roots = vec![root];
                if let Some((second, third)) =
                    quadratic_roots::<B>(candidate, &deflated_linear, &deflated_constant)?
                {
                    roots.push(second);
                    roots.push(third);
                }
                roots
            }
            3 => {
                let inverse = factor[2].modular_inverse(candidate)?;
                let monic_linear = modular_mul::<B>(&factor[1], &inverse, candidate)?;
                let monic_constant = modular_mul::<B>(&factor[0], &inverse, candidate)?;
                let Some((first, second)) =
                    quadratic_roots::<B>(candidate, &monic_linear, &monic_constant)?
                else {
                    continue;
                };
                // The root sum is −quadratic, which yields the third root.
                let sum = ((first.clone() + &second)? % candidate)?;
                let third = modular_sub::<B>(
                    &negate_mod::<B>(&polynomial.quadratic, candidate)?,
                    &sum,
                    candidate,
                )?;
                vec![first, second, third]
            }
            _ => continue,
        };
        let mut roots = Vec::with_capacity(candidates.len());
        for root in candidates {
            if polynomial.evaluate(candidate, &root)?.is_zero() {
                roots.push(root);
            }
        }
        if !roots.is_empty() {
            return Ok(Some(roots));
        }
    }
    Ok(None)
}

fn j_invariants<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    candidate: &B,
    polynomial: ClassPolynomial,
    rng: &mut R,
) -> Result<Option<Vec<B>>> {
    match polynomial {
        ClassPolynomial::Linear(root) => Ok(Some(vec![modular_signed::<B>(root, candidate)?])),
        ClassPolynomial::Quadratic { constant, linear } => {
            let linear_mod = modular_signed::<B>(linear, candidate)?;
            let constant_mod = modular_signed::<B>(constant, candidate)?;
            Ok(quadratic_roots::<B>(candidate, &linear_mod, &constant_mod)?
                .map(|(first, second)| vec![first, second]))
        }
        ClassPolynomial::Cubic {
            constant,
            linear,
            quadratic,
        } => {
            let polynomial = CubicPolynomial::<B>::new(
                candidate,
                modular_signed::<B>(quadratic, candidate)?,
                modular_signed::<B>(linear, candidate)?,
                modular_signed::<B>(constant, candidate)?,
            )?;
            cubic_roots::<B, R>(candidate, &polynomial, rng)
        }
    }
}

fn cornacchia<B: ArithmeticBackend>(candidate: &B, discriminant: i16) -> Result<Option<(B, B)>> {
    let absolute = u64::from(discriminant.unsigned_abs());
    let residue = modular_signed::<B>(i128::from(discriminant), candidate)?;
    let Some(mut root) = modular_sqrt::<B>(&residue, candidate)? else {
        return Ok(None);
    };
    let expected_odd = absolute & 1 == 1;
    if root.bit(0) != expected_odd {
        root = (candidate.clone() - &root)?;
    }
    let mut previous = (candidate.clone() << 1)?;
    let mut current = root;
    let four_candidate = (candidate.clone() << 2)?;
    let limit = four_candidate.sqrt()?;
    while current > limit {
        let remainder = (previous % &current)?;
        previous = current;
        current = remainder;
        if current.is_zero() {
            return Ok(None);
        }
    }
    let square = (current.clone() * &current)?;
    if square > four_candidate {
        return Ok(None);
    }
    let remainder = (four_candidate - &square)?;
    let absolute = from_u64::<B>(absolute)?;
    if !(remainder.clone() % &absolute)?.is_zero() {
        return Ok(None);
    }
    let v_squared = (remainder / &absolute)?;
    let Some(v) = is_square::<B>(&v_squared)? else {
        return Ok(None);
    };
    Ok(Some((current, v)))
}

fn ceil_fourth_root<B: ArithmeticBackend>(value: &B) -> Result<B> {
    let mut root = value.sqrt()?.sqrt()?;
    let square = (root.clone() * &root)?;
    if (square.clone() * &square)? < *value {
        root = add_u64::<B>(&root, 1)?;
    }
    Ok(root)
}

#[derive(Clone)]
struct AffineCurve<B: ArithmeticBackend> {
    a: B,
    b: B,
}

enum AffinePoint<B: ArithmeticBackend> {
    Infinity,
    Finite { x: B, y: B },
}

impl<B: ArithmeticBackend> Clone for AffinePoint<B> {
    fn clone(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Finite { x, y } => Self::Finite {
                x: x.clone(),
                y: y.clone(),
            },
        }
    }
}

impl<B: ArithmeticBackend> PartialEq for AffinePoint<B> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Infinity, Self::Infinity) => true,
            (Self::Finite { x: x1, y: y1 }, Self::Finite { x: x2, y: y2 }) => x1 == x2 && y1 == y2,
            _ => false,
        }
    }
}

impl<B: ArithmeticBackend> Eq for AffinePoint<B> {}

fn point_add<B: ArithmeticBackend>(
    curve: &AffineCurve<B>,
    modulus: &B,
    left: &AffinePoint<B>,
    right: &AffinePoint<B>,
) -> Result<AffinePoint<B>> {
    let (x1, y1) = match left {
        AffinePoint::Infinity => return Ok(right.clone()),
        AffinePoint::Finite { x, y } => (x, y),
    };
    let (x2, y2) = match right {
        AffinePoint::Infinity => return Ok(left.clone()),
        AffinePoint::Finite { x, y } => (x, y),
    };

    let slope = if x1 == x2 {
        if ((y1.clone() + y2)? % modulus)?.is_zero() {
            return Ok(AffinePoint::Infinity);
        }
        if y1 != y2 {
            return Err(Error::Composite);
        }
        let denominator = (mul_u64::<B>(y1, 2)? % modulus)?;
        let x_squared = (x1.clone() * x1)?;
        let numerator = (mul_u64::<B>(&x_squared, 3)? + &curve.a)?;
        modular_mul::<B>(&numerator, &denominator.modular_inverse(modulus)?, modulus)?
    } else {
        let numerator = modular_sub::<B>(y2, y1, modulus)?;
        let denominator = modular_sub::<B>(x2, x1, modulus)?;
        modular_mul::<B>(&numerator, &denominator.modular_inverse(modulus)?, modulus)?
    };
    let slope_squared = modular_mul::<B>(&slope, &slope, modulus)?;
    let x3 = modular_sub::<B>(&modular_sub::<B>(&slope_squared, x1, modulus)?, x2, modulus)?;
    let y3 = modular_sub::<B>(
        &modular_mul::<B>(&slope, &modular_sub::<B>(x1, &x3, modulus)?, modulus)?,
        y1,
        modulus,
    )?;
    Ok(AffinePoint::Finite { x: x3, y: y3 })
}

fn scalar_mul<B: ArithmeticBackend>(
    curve: &AffineCurve<B>,
    modulus: &B,
    scalar: &B,
    point: &AffinePoint<B>,
) -> Result<AffinePoint<B>> {
    let mut output = AffinePoint::Infinity;
    let mut addend = point.clone();
    let bits = scalar.bit_length();
    for index in 0..bits {
        if scalar.bit(index) {
            output = point_add::<B>(curve, modulus, &output, &addend)?;
        }
        if index + 1 < bits {
            addend = point_add::<B>(curve, modulus, &addend, &addend)?;
        }
    }
    Ok(output)
}

/// Bounds for heap-backed certificate construction.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProverOptions {
    /// The trial-division bound used to find the smooth part of a curve order.
    pub trial_division_limit: u32,
    /// The number of random x-coordinates tried for each candidate curve.
    pub point_attempts: u32,
    /// The maximum number of ECPP reductions in one certificate.
    pub max_depth: usize,
    /// The largest CM class number searched, between one and three. Larger
    /// class numbers widen the search at the cost of polynomial root
    /// splitting work per level.
    pub max_class_number: u8,
    /// The number of ECM escalation rounds tried, up to three, when Pollard
    /// rho fails to split a curve-order cofactor of at least 384 bits. Each
    /// round raises the stage-one bound and curve count, extending the
    /// reachable factor range at growing cost. Zero disables ECM.
    pub ecm_rounds: u8,
}

impl Default for ProverOptions {
    fn default() -> Self {
        Self {
            trial_division_limit: 10_000,
            point_attempts: 128,
            max_depth: 64,
            max_class_number: 3,
            ecm_rounds: 1,
        }
    }
}

struct WorkingStep<B: ArithmeticBackend> {
    n: B,
    curve: AffineCurve<B>,
    point: AffinePoint<B>,
    cofactor: B,
    q: B,
}

pub(crate) fn prove<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    candidate: &[u8],
    rng: &mut R,
    options: ProverOptions,
) -> Result<PrimalityProof> {
    let candidate = B::from_be_bytes(candidate)?;
    if cmp_u64::<B>(&candidate, 2)? == Ordering::Less {
        return Err(Error::Composite);
    }
    let primes = small_primes(options.trial_division_limit.max(53));
    if !is_probable_prime::<B>(&candidate, &primes)? {
        return Err(Error::Composite);
    }

    let mut nodes = Vec::new();
    let mut current = candidate;
    for _ in 0..options.max_depth {
        if let Some(small) = to_u64::<B>(&current) {
            if !is_prime_u64(small) {
                return Err(Error::Composite);
            }
            nodes.push(ProofNode::SmallPrime(small));
            return Ok(PrimalityProof { nodes });
        }
        let step = find_step::<B, R>(&current, rng, &primes, options)?;
        current = step.q.clone();
        nodes.push(ProofNode::EllipticCurve(encode_step::<B>(&step)?));
    }
    Err(Error::SearchExhausted {
        candidate: crate::Natural::from_be_bytes(&current.to_be_bytes()),
    })
}

fn encode_step<B: ArithmeticBackend>(step: &WorkingStep<B>) -> Result<crate::EcppStep> {
    let point = match &step.point {
        AffinePoint::Finite { x, y } => crate::Point {
            x: crate::Natural::from_be_bytes(&x.to_be_bytes()),
            y: crate::Natural::from_be_bytes(&y.to_be_bytes()),
        },
        AffinePoint::Infinity => {
            return Err(Error::Arithmetic(
                "prover produced an infinite certificate point",
            ));
        }
    };
    Ok(crate::EcppStep {
        n: crate::Natural::from_be_bytes(&step.n.to_be_bytes()),
        curve: crate::Curve {
            a: crate::Natural::from_be_bytes(&step.curve.a.to_be_bytes()),
            b: crate::Natural::from_be_bytes(&step.curve.b.to_be_bytes()),
        },
        point,
        cofactor: crate::Natural::from_be_bytes(&step.cofactor.to_be_bytes()),
        q: crate::Natural::from_be_bytes(&step.q.to_be_bytes()),
    })
}

struct StepSearch<'a, B: ArithmeticBackend, R: ?Sized> {
    candidate: &'a B,
    rng: &'a mut R,
    primes: &'a [u32],
    point_attempts: u32,
    ecm_rounds: u8,
}

fn find_step<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    candidate: &B,
    rng: &mut R,
    primes: &[u32],
    options: ProverOptions,
) -> Result<WorkingStep<B>> {
    let mut search = StepSearch::<B, R> {
        candidate,
        rng,
        primes,
        point_attempts: options.point_attempts,
        ecm_rounds: options.ecm_rounds,
    };
    let order = add_u64::<B>(candidate, 1)?;
    // A value is 3 modulo 4 exactly when its low two bits are set.
    if candidate.bit(0) && candidate.bit(1) {
        let curve = AffineCurve::<B> {
            a: one::<B>()?,
            b: zero::<B>()?,
        };
        if let Some(step) = try_order::<B, R>(&mut search, &order, &curve, None)? {
            return Ok(step);
        }
    }
    let modulo_three = (candidate.clone() % &from_u64::<B>(3)?)?;
    if cmp_u64::<B>(&modulo_three, 2)? == Ordering::Equal {
        let curve = AffineCurve::<B> {
            a: zero::<B>()?,
            b: one::<B>()?,
        };
        if let Some(step) = try_order::<B, R>(&mut search, &order, &curve, None)? {
            return Ok(step);
        }
    }

    for discriminant in DISCRIMINANTS {
        if discriminant.polynomial.class_number() > options.max_class_number {
            continue;
        }
        let Some((trace, _)) = cornacchia::<B>(candidate, discriminant.value)? else {
            continue;
        };
        let Some(invariants) =
            j_invariants::<B, R>(candidate, discriminant.polynomial, &mut *search.rng)?
        else {
            continue;
        };
        for invariant in invariants {
            let base = curve_from_j::<B>(candidate, &invariant)?;
            let twist = quadratic_twist::<B>(candidate, &base)?;
            let lower = (order.clone() - &trace)?;
            let upper = (order.clone() + &trace)?;
            for curve_order in [lower, upper] {
                if let Some(step) =
                    try_order::<B, R>(&mut search, &curve_order, &base, Some(&twist))?
                {
                    return Ok(step);
                }
            }
        }
    }
    Err(Error::SearchExhausted {
        candidate: crate::Natural::from_be_bytes(&candidate.to_be_bytes()),
    })
}

fn try_order<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    search: &mut StepSearch<'_, B, R>,
    order: &B,
    curve: &AffineCurve<B>,
    twist: Option<&AffineCurve<B>>,
) -> Result<Option<WorkingStep<B>>> {
    let Some(q) = split_order::<B, R>(
        search.candidate,
        order,
        search.primes,
        search.ecm_rounds,
        search.rng,
    )?
    else {
        return Ok(None);
    };
    let cofactor = (order.clone() / &q)?;
    for candidate_curve in core::iter::once(curve).chain(twist) {
        if let Some(point) = find_point_of_order::<B, R>(
            search.candidate,
            candidate_curve,
            order,
            &q,
            search.rng,
            search.point_attempts,
        )? {
            return Ok(Some(WorkingStep {
                n: search.candidate.clone(),
                curve: AffineCurve {
                    a: candidate_curve.a.clone(),
                    b: candidate_curve.b.clone(),
                },
                point,
                cofactor,
                q,
            }));
        }
    }
    Ok(None)
}

fn split_order<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    candidate: &B,
    order: &B,
    primes: &[u32],
    ecm_rounds: u8,
    rng: &mut R,
) -> Result<Option<B>> {
    let mut remaining = order.clone();
    for &prime in primes {
        let prime = from_u64::<B>(u64::from(prime))?;
        while (remaining.clone() % &prime)?.is_zero() {
            remaining = (remaining / &prime)?;
        }
    }
    let fourth_root = ceil_fourth_root::<B>(candidate)?;
    let root_plus_one = add_u64::<B>(&fourth_root, 1)?;
    let bound = (root_plus_one.clone() * &root_plus_one)?;
    find_large_prime_factor::<B, R>(&remaining, candidate, &bound, primes, ecm_rounds, rng, 0)
}

/// The smallest cofactor width where ECM escalation is attempted after
/// Pollard rho fails; below it, rho already covers the reachable factors.
const ECM_MINIMUM_BITS: usize = 384;

fn find_large_prime_factor<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    value: &B,
    candidate: &B,
    bound: &B,
    primes: &[u32],
    ecm_rounds: u8,
    rng: &mut R,
    depth: usize,
) -> Result<Option<B>> {
    if value <= bound || depth > 32 {
        return Ok(None);
    }
    if is_probable_prime::<B>(value, primes)? {
        return Ok((value < candidate).then(|| value.clone()));
    }
    let factor = match pollard_rho::<B, R>(value, rng)? {
        Some(factor) => factor,
        None if ecm_rounds > 0 && value.bit_length() >= ECM_MINIMUM_BITS => {
            let Some(factor) = ecm::<B, R>(value, ecm_rounds, rng)? else {
                return Ok(None);
            };
            factor
        }
        None => return Ok(None),
    };
    let other = (value.clone() / &factor)?;
    let (first, second) = if factor >= other {
        (factor, other)
    } else {
        (other, factor)
    };
    if let Some(found) = find_large_prime_factor::<B, R>(
        &first,
        candidate,
        bound,
        primes,
        ecm_rounds,
        rng,
        depth + 1,
    )? {
        return Ok(Some(found));
    }
    find_large_prime_factor::<B, R>(
        &second,
        candidate,
        bound,
        primes,
        ecm_rounds,
        rng,
        depth + 1,
    )
}

/// A projective x-only point on a Montgomery curve modulo the composite
/// being factored.
struct MontgomeryPoint<B: ArithmeticBackend> {
    x: B,
    z: B,
}

impl<B: ArithmeticBackend> Clone for MontgomeryPoint<B> {
    fn clone(&self) -> Self {
        Self {
            x: self.x.clone(),
            z: self.z.clone(),
        }
    }
}

/// x-only Montgomery curve arithmetic modulo the composite being factored,
/// parameterized by `a24 = (A + 2) / 4`.
struct EcmCurve<'a, B: ArithmeticBackend> {
    value: &'a B,
    a24: B,
}

impl<B: ArithmeticBackend> EcmCurve<'_, B> {
    fn double(&self, point: &MontgomeryPoint<B>) -> Result<MontgomeryPoint<B>> {
        let sum = ((point.x.clone() + &point.z)? % self.value)?;
        let sum_squared = modular_mul::<B>(&sum, &sum, self.value)?;
        let difference = modular_sub::<B>(&point.x, &point.z, self.value)?;
        let difference_squared = modular_mul::<B>(&difference, &difference, self.value)?;
        let x = modular_mul::<B>(&sum_squared, &difference_squared, self.value)?;
        let four_xz = modular_sub::<B>(&sum_squared, &difference_squared, self.value)?;
        let scaled = modular_mul::<B>(&self.a24, &four_xz, self.value)?;
        let z = modular_mul::<B>(
            &four_xz,
            &((difference_squared + &scaled)? % self.value)?,
            self.value,
        )?;
        Ok(MontgomeryPoint { x, z })
    }

    fn add(
        &self,
        left: &MontgomeryPoint<B>,
        right: &MontgomeryPoint<B>,
        difference: &MontgomeryPoint<B>,
    ) -> Result<MontgomeryPoint<B>> {
        let left_minus = modular_sub::<B>(&left.x, &left.z, self.value)?;
        let left_plus = ((left.x.clone() + &left.z)? % self.value)?;
        let right_minus = modular_sub::<B>(&right.x, &right.z, self.value)?;
        let right_plus = ((right.x.clone() + &right.z)? % self.value)?;
        let cross_one = modular_mul::<B>(&left_minus, &right_plus, self.value)?;
        let cross_two = modular_mul::<B>(&left_plus, &right_minus, self.value)?;
        let sum = ((cross_one.clone() + &cross_two)? % self.value)?;
        let difference_of_crosses = modular_sub::<B>(&cross_one, &cross_two, self.value)?;
        let x = modular_mul::<B>(
            &difference.z,
            &modular_mul::<B>(&sum, &sum, self.value)?,
            self.value,
        )?;
        let z = modular_mul::<B>(
            &difference.x,
            &modular_mul::<B>(&difference_of_crosses, &difference_of_crosses, self.value)?,
            self.value,
        )?;
        Ok(MontgomeryPoint { x, z })
    }

    fn ladder(&self, scalar: u32, point: &MontgomeryPoint<B>) -> Result<MontgomeryPoint<B>> {
        if scalar == 0 {
            return Ok(MontgomeryPoint {
                x: one::<B>()?,
                z: zero::<B>()?,
            });
        }
        if scalar == 1 {
            return Ok(point.clone());
        }
        let mut lower = point.clone();
        let mut upper = self.double(point)?;
        let bits = 32 - scalar.leading_zeros();
        for index in (0..bits - 1).rev() {
            if scalar & (1 << index) != 0 {
                lower = self.add(&lower, &upper, point)?;
                upper = self.double(&upper)?;
            } else {
                upper = self.add(&lower, &upper, point)?;
                lower = self.double(&lower)?;
            }
        }
        Ok(lower)
    }
}

/// ECM escalation rounds as `(stage-one bound, curve count)`. Pollard rho
/// already covers factors reachable below the first round's bound.
const ECM_ROUNDS: [(u32, u32); 3] = [(11_000, 20), (50_000, 40), (250_000, 80)];

/// Runs one Suyama-parameterized ECM curve and returns a nontrivial factor
/// when stage one finds one.
fn ecm_curve_attempt<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    value: &B,
    primes: &[u32],
    bound1: u32,
    rng: &mut R,
) -> Result<Option<B>> {
    let sigma = random_below::<B, R>(value, rng)?;
    if cmp_u64::<B>(&sigma, 6)? == Ordering::Less {
        return Ok(None);
    }
    let sigma_squared = modular_mul::<B>(&sigma, &sigma, value)?;
    let u = modular_sub::<B>(&sigma_squared, &from_u64::<B>(5)?, value)?;
    let v = (mul_u64::<B>(&sigma, 4)? % value)?;
    if u.is_zero() || v.is_zero() {
        return Ok(None);
    }
    let u_squared = modular_mul::<B>(&u, &u, value)?;
    let u_cubed = modular_mul::<B>(&u_squared, &u, value)?;
    let v_squared = modular_mul::<B>(&v, &v, value)?;
    let v_cubed = modular_mul::<B>(&v_squared, &v, value)?;

    // a24 = (v − u)³ (3u + v) / (16 u³ v) for the curve constant.
    let difference = modular_sub::<B>(&v, &u, value)?;
    let difference_squared = modular_mul::<B>(&difference, &difference, value)?;
    let difference_cubed = modular_mul::<B>(&difference_squared, &difference, value)?;
    let three_u_plus_v = ((mul_u64::<B>(&u, 3)? + &v)? % value)?;
    let numerator = modular_mul::<B>(&difference_cubed, &three_u_plus_v, value)?;
    let denominator = (mul_u64::<B>(&modular_mul::<B>(&u_cubed, &v, value)?, 16)? % value)?;
    let shared = denominator.gcd(value)?;
    if !shared.is_one() {
        // A degenerate parameterization either exposes a factor directly or
        // collapses entirely; only the former is usable.
        return Ok((&shared != value).then_some(shared));
    }
    let a24 = modular_mul::<B>(&numerator, &denominator.modular_inverse(value)?, value)?;

    let curve = EcmCurve { value, a24 };
    let mut point = MontgomeryPoint {
        x: u_cubed,
        z: v_cubed,
    };
    let limit = u64::from(bound1);
    for &prime in primes {
        let mut power = u64::from(prime);
        loop {
            point = curve.ladder(prime, &point)?;
            power = power.saturating_mul(u64::from(prime));
            if power > limit {
                break;
            }
        }
    }
    let shared = point.z.gcd(value)?;
    Ok((!shared.is_one() && &shared != value).then_some(shared))
}

/// Lenstra's elliptic curve method, tried after Pollard rho fails on a
/// large cofactor. Rounds escalate the stage-one bound per `ECM_ROUNDS`.
fn ecm<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    value: &B,
    rounds: u8,
    rng: &mut R,
) -> Result<Option<B>> {
    for (bound1, curves) in ECM_ROUNDS.iter().take(usize::from(rounds)) {
        let primes = small_primes(*bound1);
        for _ in 0..*curves {
            if let Some(factor) = ecm_curve_attempt::<B, R>(value, &primes, *bound1, rng)? {
                return Ok(Some(factor));
            }
        }
    }
    Ok(None)
}

fn pollard_rho<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    value: &B,
    rng: &mut R,
) -> Result<Option<B>> {
    if value.is_even() {
        return Ok(Some(from_u64::<B>(2)?));
    }
    let one = one::<B>()?;
    const BATCH: usize = 128;
    const MAX_ITERATIONS: usize = 1_000_000;
    for _ in 0..16 {
        let mut y = random_below::<B, R>(value, rng)?;
        let constant = add_u64::<B>(&random_below::<B, R>(value, rng)?, 1)?;
        let mut divisor = one.clone();
        let mut power = 1usize;
        let mut iterations = 0usize;
        let mut x = y.clone();
        let mut saved_y = y.clone();
        while divisor.is_one() && iterations < MAX_ITERATIONS {
            x = y.clone();
            for _ in 0..power {
                y = (((y.clone() * &y)? + &constant)? % value)?;
            }
            iterations += power;
            let mut offset = 0usize;
            while offset < power && divisor.is_one() {
                saved_y = y.clone();
                let count = BATCH.min(power - offset);
                let mut product = one.clone();
                for _ in 0..count {
                    y = (((y.clone() * &y)? + &constant)? % value)?;
                    let difference = if x >= y {
                        (x.clone() - &y)?
                    } else {
                        (y.clone() - &x)?
                    };
                    product = ((product * &difference)? % value)?;
                }
                divisor = product.gcd(value)?;
                offset += count;
                iterations += count;
            }
            power = power.saturating_mul(2);
        }
        if divisor == *value {
            divisor = one.clone();
            while divisor.is_one() && iterations < MAX_ITERATIONS * 2 {
                saved_y = (((saved_y.clone() * &saved_y)? + &constant)? % value)?;
                let difference = if x >= saved_y {
                    (x.clone() - &saved_y)?
                } else {
                    (saved_y.clone() - &x)?
                };
                divisor = difference.gcd(value)?;
                iterations += 1;
            }
        }
        if !divisor.is_one() && divisor != *value {
            return Ok(Some(divisor));
        }
    }
    Ok(None)
}

fn curve_from_j<B: ArithmeticBackend>(candidate: &B, invariant: &B) -> Result<AffineCurve<B>> {
    let denominator = modular_sub::<B>(&from_u64::<B>(1728)?, invariant, candidate)?;
    let inverse = denominator.modular_inverse(candidate)?;
    let k = modular_mul::<B>(invariant, &inverse, candidate)?;
    Ok(AffineCurve {
        a: (mul_u64::<B>(&k, 3)? % candidate)?,
        b: (mul_u64::<B>(&k, 2)? % candidate)?,
    })
}

fn quadratic_twist<B: ArithmeticBackend>(
    candidate: &B,
    curve: &AffineCurve<B>,
) -> Result<AffineCurve<B>> {
    let mut non_residue = from_u64::<B>(2)?;
    while non_residue.jacobi(candidate)? != -1 {
        non_residue = add_u64::<B>(&non_residue, 1)?;
        if &non_residue >= candidate {
            return Err(Error::Composite);
        }
    }
    let square = modular_mul::<B>(&non_residue, &non_residue, candidate)?;
    Ok(AffineCurve {
        a: modular_mul::<B>(&curve.a, &square, candidate)?,
        b: modular_mul::<B>(
            &modular_mul::<B>(&curve.b, &square, candidate)?,
            &non_residue,
            candidate,
        )?,
    })
}

fn find_point_of_order<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    modulus: &B,
    curve: &AffineCurve<B>,
    order: &B,
    q: &B,
    rng: &mut R,
    attempts: u32,
) -> Result<Option<AffinePoint<B>>> {
    let cofactor = (order.clone() / q)?;
    for _ in 0..attempts {
        let x = random_below::<B, R>(modulus, rng)?;
        let x_squared = modular_mul::<B>(&x, &x, modulus)?;
        let x_cubed = modular_mul::<B>(&x_squared, &x, modulus)?;
        let ax = modular_mul::<B>(&curve.a, &x, modulus)?;
        let rhs = (((x_cubed + &ax)? + &curve.b)? % modulus)?;
        let Some(y) = modular_sqrt::<B>(&rhs, modulus)? else {
            continue;
        };
        let point = AffinePoint::Finite { x, y };
        let q_point = scalar_mul::<B>(curve, modulus, &cofactor, &point)?;
        if q_point == AffinePoint::Infinity {
            continue;
        }
        if scalar_mul::<B>(curve, modulus, q, &q_point)? == AffinePoint::Infinity {
            return Ok(Some(point));
        }
    }
    Ok(None)
}

fn random_below<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    modulus: &B,
    rng: &mut R,
) -> Result<B> {
    let bits = modulus.bit_length();
    let byte_length = bits.div_ceil(8);
    let excess = byte_length * 8 - bits;
    loop {
        let mut bytes = vec![0u8; byte_length];
        rng.fill_bytes(&mut bytes);
        bytes[0] &= u8::MAX >> excess;
        let value = B::from_be_bytes(&bytes)?;
        if &value < modulus {
            return Ok(value);
        }
    }
}

pub(crate) fn verify_proof<B: ArithmeticBackend>(proof: &PrimalityProof) -> Result<()> {
    if proof.nodes.is_empty() {
        return Err(Error::InvalidProof("certificate is empty"));
    }
    let mut expected: Option<B> = None;
    for (index, node) in proof.nodes.iter().enumerate() {
        match node {
            ProofNode::SmallPrime(prime) => {
                if index + 1 != proof.nodes.len() {
                    return Err(Error::InvalidProof("small-prime node must be last"));
                }
                let base_case = from_u64::<B>(*prime)?;
                if expected.as_ref().is_some_and(|value| value != &base_case) {
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
                verify_step::<B>(step, expected.as_ref())?;
                expected = Some(B::from_be_bytes(step.q.as_be_bytes())?);
            }
        }
    }
    Ok(())
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

fn verify_step<B: ArithmeticBackend>(step: &crate::EcppStep, expected: Option<&B>) -> Result<()> {
    let n = B::from_be_bytes(step.n.as_be_bytes())?;
    let q = B::from_be_bytes(step.q.as_be_bytes())?;
    let cofactor = B::from_be_bytes(step.cofactor.as_be_bytes())?;
    if expected.is_some_and(|value| value != &n) {
        return Err(Error::InvalidProof("certificate chain is disconnected"));
    }
    if n.is_even()
        || cmp_u64::<B>(&n, 3)? == Ordering::Less
        || q >= n
        || cmp_u64::<B>(&q, 2)? == Ordering::Less
    {
        return Err(Error::InvalidProof("invalid ECPP step integers"));
    }
    let fourth_root = ceil_fourth_root::<B>(&n)?;
    let root_plus_one = add_u64::<B>(&fourth_root, 1)?;
    if q <= (root_plus_one.clone() * &root_plus_one)? {
        return Err(Error::InvalidProof(
            "q is below the elliptic Pocklington bound",
        ));
    }

    let curve = AffineCurve::<B> {
        a: B::from_be_bytes(step.curve.a.as_be_bytes())?,
        b: B::from_be_bytes(step.curve.b.as_be_bytes())?,
    };
    let point = AffinePoint::<B>::Finite {
        x: B::from_be_bytes(step.point.x.as_be_bytes())?,
        y: B::from_be_bytes(step.point.y.as_be_bytes())?,
    };
    if curve.a >= n || curve.b >= n {
        return Err(Error::InvalidProof("curve coefficients are not reduced"));
    }
    let (x, y) = match &point {
        AffinePoint::Finite { x, y } if x < &n && y < &n => (x, y),
        _ => return Err(Error::InvalidProof("point coordinates are not reduced")),
    };
    let a_squared = modular_mul::<B>(&curve.a, &curve.a, &n)?;
    let a_cubed = modular_mul::<B>(&a_squared, &curve.a, &n)?;
    let b_squared = modular_mul::<B>(&curve.b, &curve.b, &n)?;
    let discriminant = ((mul_u64::<B>(&a_cubed, 4)? + &mul_u64::<B>(&b_squared, 27)?)? % &n)?;
    if !discriminant.gcd(&n)?.is_one() {
        return Err(Error::InvalidProof(
            "curve is singular modulo a divisor of n",
        ));
    }
    let x_squared = modular_mul::<B>(x, x, &n)?;
    let x_cubed = modular_mul::<B>(&x_squared, x, &n)?;
    let ax = modular_mul::<B>(&curve.a, x, &n)?;
    let rhs = (((x_cubed + &ax)? + &curve.b)? % &n)?;
    if modular_mul::<B>(y, y, &n)? != rhs {
        return Err(Error::InvalidProof("certificate point is not on the curve"));
    }

    let q_point = scalar_mul::<B>(&curve, &n, &cofactor, &point)
        .map_err(|_| Error::InvalidProof("point multiplication is not defined modulo n"))?;
    if q_point == AffinePoint::Infinity {
        return Err(Error::InvalidProof(
            "cofactor annihilates the certificate point",
        ));
    }
    let result = scalar_mul::<B>(&curve, &n, &q, &q_point)
        .map_err(|_| Error::InvalidProof("point multiplication is not defined modulo n"))?;
    if result != AffinePoint::Infinity {
        return Err(Error::InvalidProof(
            "curve order does not annihilate the point",
        ));
    }
    Ok(())
}

#[cfg(all(test, feature = "num-bigint"))]
mod tests {
    use super::*;
    use crate::arithmetic::NumBigint;
    use rand::{SeedableRng, rngs::StdRng};

    #[test]
    fn default_options_search_the_full_table() {
        assert_eq!(ProverOptions::default().max_class_number, 3);
        assert_eq!(ProverOptions::default().ecm_rounds, 1);
    }

    #[test]
    fn ecm_splits_a_semiprime_beyond_trial_range() {
        let mut rng = StdRng::seed_from_u64(9);
        let first = from_u64::<NumBigint>(1_000_000_007).unwrap();
        let second = from_u64::<NumBigint>(1_000_000_009).unwrap();
        let value = (first.clone() * &second).unwrap();
        let factor = ecm::<NumBigint, _>(&value, 2, &mut rng).unwrap().unwrap();
        assert!(factor == first || factor == second);
    }

    #[test]
    fn cubic_class_polynomial_splits_for_represented_primes() {
        let cubic_entry = DISCRIMINANTS
            .iter()
            .find(|discriminant| discriminant.value == -23)
            .expect("table contains -23");
        let ClassPolynomial::Cubic {
            constant,
            linear,
            quadratic,
        } = cubic_entry.polynomial
        else {
            panic!("-23 should be cubic");
        };

        let mut rng = StdRng::seed_from_u64(5);
        // 4·59 = 12² + 23·2² and 4·101 = 6² + 23·4², so H₋₂₃ splits
        // completely modulo both primes.
        for prime in [59u64, 101] {
            let candidate = from_u64::<NumBigint>(prime).unwrap();
            let polynomial = CubicPolynomial::<NumBigint>::new(
                &candidate,
                modular_signed::<NumBigint>(quadratic, &candidate).unwrap(),
                modular_signed::<NumBigint>(linear, &candidate).unwrap(),
                modular_signed::<NumBigint>(constant, &candidate).unwrap(),
            )
            .unwrap();
            let roots = cubic_roots::<NumBigint, _>(&candidate, &polynomial, &mut rng)
                .unwrap()
                .unwrap();
            assert_eq!(roots.len(), 3, "prime {prime} should split completely");
            for root in &roots {
                assert!(polynomial.evaluate(&candidate, root).unwrap().is_zero());
            }
        }
    }
}
