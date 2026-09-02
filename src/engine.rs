use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

use rand_core::CryptoRng;

use crate::arithmetic::{
    ArithmeticBackend, bit, bit_length, cmp_u64, from_u64, from_u128, is_even, to_u64,
};
use crate::cm::{ClassPolynomial, DISCRIMINANTS};
use crate::{Error, PrimalityProof, ProofNode, Result};

fn add<B: ArithmeticBackend>(left: &B, right: &B) -> Result<B> {
    left.clone() + right
}

fn sub<B: ArithmeticBackend>(left: &B, right: &B) -> Result<B> {
    left.clone() - right
}

fn mul<B: ArithmeticBackend>(left: &B, right: &B) -> Result<B> {
    left.clone() * right
}

fn div<B: ArithmeticBackend>(left: &B, right: &B) -> Result<B> {
    left.clone() / right
}

fn rem<B: ArithmeticBackend>(left: &B, right: &B) -> Result<B> {
    left.clone() % right
}

fn shl<B: ArithmeticBackend>(value: &B, bits: usize) -> Result<B> {
    value.clone() << bits
}

fn shr<B: ArithmeticBackend>(value: &B, bits: usize) -> Result<B> {
    value.clone() >> bits
}

fn zero<B: ArithmeticBackend>() -> Result<B> {
    from_u64::<B>(0)
}

fn one<B: ArithmeticBackend>() -> Result<B> {
    from_u64::<B>(1)
}

fn is_zero<B: ArithmeticBackend>(value: &B) -> Result<bool> {
    Ok(value == &zero::<B>()?)
}

fn add_u64<B: ArithmeticBackend>(value: &B, other: u64) -> Result<B> {
    add::<B>(value, &from_u64::<B>(other)?)
}

fn mul_u64<B: ArithmeticBackend>(value: &B, other: u64) -> Result<B> {
    mul::<B>(value, &from_u64::<B>(other)?)
}

fn modular_sub<B: ArithmeticBackend>(left: &B, right: &B, modulus: &B) -> Result<B> {
    if left >= right {
        rem::<B>(&sub::<B>(left, right)?, modulus)
    } else {
        let difference = rem::<B>(&sub::<B>(right, left)?, modulus)?;
        if is_zero::<B>(&difference)? {
            zero::<B>()
        } else {
            sub::<B>(modulus, &difference)
        }
    }
}

fn modular_mul<B: ArithmeticBackend>(left: &B, right: &B, modulus: &B) -> Result<B> {
    rem::<B>(&mul::<B>(left, right)?, modulus)
}

fn modular_pow<B: ArithmeticBackend>(base: &B, exponent: &B, modulus: &B) -> Result<B> {
    let mut output = one::<B>()?;
    let mut base = rem::<B>(base, modulus)?;
    let mut index = 0;
    while index < crate::arithmetic::bit_length::<B>(exponent) {
        if bit::<B>(exponent, index) {
            output = modular_mul::<B>(&output, &base, modulus)?;
        }
        base = modular_mul::<B>(&base, &base, modulus)?;
        index += 1;
    }
    Ok(output)
}

fn modular_signed<B: ArithmeticBackend>(value: i128, modulus: &B) -> Result<B> {
    let magnitude = rem::<B>(&from_u128::<B>(value.unsigned_abs())?, modulus)?;
    if value >= 0 || is_zero::<B>(&magnitude)? {
        Ok(magnitude)
    } else {
        sub::<B>(modulus, &magnitude)
    }
}

fn gcd<B: ArithmeticBackend>(left: &B, right: &B) -> Result<B> {
    let mut left = left.clone();
    let mut right = right.clone();
    while !is_zero::<B>(&right)? {
        let remainder = rem::<B>(&left, &right)?;
        left = right;
        right = remainder;
    }
    Ok(left)
}

fn modular_inverse<B: ArithmeticBackend>(value: &B, modulus: &B) -> Result<B> {
    let mut coefficient = zero::<B>()?;
    let mut next_coefficient = one::<B>()?;
    let mut remainder = modulus.clone();
    let mut next_remainder = rem::<B>(value, modulus)?;
    while !is_zero::<B>(&next_remainder)? {
        let quotient = div::<B>(&remainder, &next_remainder)?;
        let old_coefficient = coefficient;
        coefficient = next_coefficient;
        let product = rem::<B>(&mul::<B>(&quotient, &coefficient)?, modulus)?;
        next_coefficient = modular_sub::<B>(&old_coefficient, &product, modulus)?;
        let old_remainder = remainder;
        let following_remainder = sub::<B>(&old_remainder, &mul::<B>(&quotient, &next_remainder)?)?;
        remainder = next_remainder;
        next_remainder = following_remainder;
    }
    if remainder != one::<B>()? {
        return Err(Error::Composite);
    }
    Ok(coefficient)
}

fn integer_sqrt<B: ArithmeticBackend>(value: &B) -> Result<B> {
    if cmp_u64::<B>(value, 2)? == Ordering::Less {
        return Ok(value.clone());
    }
    let shift = crate::arithmetic::bit_length::<B>(value).div_ceil(2);
    let mut estimate = shl::<B>(&one::<B>()?, shift)?;
    loop {
        let next = shr::<B>(&add::<B>(&estimate, &div::<B>(value, &estimate)?)?, 1)?;
        if next >= estimate {
            return Ok(estimate);
        }
        estimate = next;
    }
}

fn is_square<B: ArithmeticBackend>(value: &B) -> Result<Option<B>> {
    let root = integer_sqrt::<B>(value)?;
    Ok((mul::<B>(&root, &root)? == *value).then_some(root))
}

fn jacobi<B: ArithmeticBackend>(value: &B, modulus: &B) -> Result<i8> {
    if is_zero::<B>(modulus)? || is_even::<B>(modulus) {
        return Ok(0);
    }
    let mut value = rem::<B>(value, modulus)?;
    let mut modulus = modulus.clone();
    let mut result = 1i8;
    while !is_zero::<B>(&value)? {
        while is_even::<B>(&value) {
            value = shr::<B>(&value, 1)?;
            let reduced = rem::<B>(&modulus, &from_u64::<B>(8)?)?;
            let residue = to_u64::<B>(&reduced).unwrap_or(0);
            if residue == 3 || residue == 5 {
                result = -result;
            }
        }
        core::mem::swap(&mut value, &mut modulus);
        let value_mod_four = rem::<B>(&value, &from_u64::<B>(4)?)?;
        let modulus_mod_four = rem::<B>(&modulus, &from_u64::<B>(4)?)?;
        if cmp_u64::<B>(&value_mod_four, 3)? == Ordering::Equal
            && cmp_u64::<B>(&modulus_mod_four, 3)? == Ordering::Equal
        {
            result = -result;
        }
        value = rem::<B>(&value, &modulus)?;
    }
    Ok(if modulus == one::<B>()? { result } else { 0 })
}

fn modular_sqrt<B: ArithmeticBackend>(value: &B, modulus: &B) -> Result<Option<B>> {
    let value = rem::<B>(value, modulus)?;
    if is_zero::<B>(&value)? {
        return Ok(Some(zero::<B>()?));
    }
    if cmp_u64::<B>(modulus, 2)? == Ordering::Equal {
        return Ok(Some(value));
    }
    if jacobi::<B>(&value, modulus)? != 1 {
        return Ok(None);
    }
    let modulus_mod_four = rem::<B>(modulus, &from_u64::<B>(4)?)?;
    if cmp_u64::<B>(&modulus_mod_four, 3)? == Ordering::Equal {
        let exponent = shr::<B>(&add_u64::<B>(modulus, 1)?, 2)?;
        return Ok(Some(modular_pow::<B>(&value, &exponent, modulus)?));
    }

    let mut odd = sub::<B>(modulus, &one::<B>()?)?;
    let mut exponent = 0u32;
    while is_even::<B>(&odd) {
        odd = shr::<B>(&odd, 1)?;
        exponent += 1;
    }
    let mut non_residue = from_u64::<B>(2)?;
    while jacobi::<B>(&non_residue, modulus)? != -1 {
        non_residue = add_u64::<B>(&non_residue, 1)?;
        if &non_residue >= modulus {
            return Ok(None);
        }
    }
    let mut c = modular_pow::<B>(&non_residue, &odd, modulus)?;
    let x_exponent = shr::<B>(&add_u64::<B>(&odd, 1)?, 1)?;
    let mut x = modular_pow::<B>(&value, &x_exponent, modulus)?;
    let mut t = modular_pow::<B>(&value, &odd, modulus)?;
    let mut m = exponent;
    while t != one::<B>()? {
        let mut i = 1u32;
        let mut power = modular_mul::<B>(&t, &t, modulus)?;
        while power != one::<B>()? {
            power = modular_mul::<B>(&power, &power, modulus)?;
            i += 1;
            if i >= m {
                return Ok(None);
            }
        }
        let two_power = shl::<B>(&one::<B>()?, (m - i - 1) as usize)?;
        let b = modular_pow::<B>(&c, &two_power, modulus)?;
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
        if is_zero::<B>(&rem::<B>(candidate, &prime_value)?)? {
            return Ok(false);
        }
    }
    if is_even::<B>(candidate) {
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
    let minus_one = sub::<B>(candidate, &one)?;
    let mut odd = minus_one.clone();
    let mut exponent = 0u32;
    while is_even::<B>(&odd) {
        odd = shr::<B>(&odd, 1)?;
        exponent += 1;
    }
    let mut value = modular_pow::<B>(base, &odd, candidate)?;
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

fn j_invariants<B: ArithmeticBackend>(
    candidate: &B,
    polynomial: ClassPolynomial,
) -> Result<Option<Vec<B>>> {
    match polynomial {
        ClassPolynomial::Linear(root) => Ok(Some(vec![modular_signed::<B>(root, candidate)?])),
        ClassPolynomial::Quadratic { constant, linear } => {
            let linear_mod = modular_signed::<B>(linear, candidate)?;
            let constant_mod = modular_signed::<B>(constant, candidate)?;
            let linear_squared = modular_mul::<B>(&linear_mod, &linear_mod, candidate)?;
            let four_constant = rem::<B>(&mul_u64::<B>(&constant_mod, 4)?, candidate)?;
            let discriminant = modular_sub::<B>(&linear_squared, &four_constant, candidate)?;
            let Some(square_root) = modular_sqrt::<B>(&discriminant, candidate)? else {
                return Ok(None);
            };
            let inverse_two = shr::<B>(&add_u64::<B>(candidate, 1)?, 1)?;
            let minus_linear = if is_zero::<B>(&linear_mod)? {
                zero::<B>()?
            } else {
                sub::<B>(candidate, &linear_mod)?
            };
            let first = modular_mul::<B>(
                &modular_sub::<B>(&minus_linear, &square_root, candidate)?,
                &inverse_two,
                candidate,
            )?;
            let second = modular_mul::<B>(
                &add::<B>(&minus_linear, &square_root)?,
                &inverse_two,
                candidate,
            )?;
            Ok(Some(vec![first, second]))
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
    if bit::<B>(&root, 0) != expected_odd {
        root = sub::<B>(candidate, &root)?;
    }
    let mut previous = shl::<B>(candidate, 1)?;
    let mut current = root;
    let four_candidate = shl::<B>(candidate, 2)?;
    let limit = integer_sqrt::<B>(&four_candidate)?;
    while current > limit {
        let remainder = rem::<B>(&previous, &current)?;
        previous = current;
        current = remainder;
        if is_zero::<B>(&current)? {
            return Ok(None);
        }
    }
    let square = mul::<B>(&current, &current)?;
    if square > four_candidate {
        return Ok(None);
    }
    let remainder = sub::<B>(&four_candidate, &square)?;
    let absolute = from_u64::<B>(absolute)?;
    if !is_zero::<B>(&rem::<B>(&remainder, &absolute)?)? {
        return Ok(None);
    }
    let v_squared = div::<B>(&remainder, &absolute)?;
    let Some(v) = is_square::<B>(&v_squared)? else {
        return Ok(None);
    };
    Ok(Some((current, v)))
}

fn ceil_fourth_root<B: ArithmeticBackend>(value: &B) -> Result<B> {
    let mut root = integer_sqrt::<B>(&integer_sqrt::<B>(value)?)?;
    let square = mul::<B>(&root, &root)?;
    if mul::<B>(&square, &square)? < *value {
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
        if is_zero::<B>(&rem::<B>(&add::<B>(y1, y2)?, modulus)?)? {
            return Ok(AffinePoint::Infinity);
        }
        if y1 != y2 {
            return Err(Error::Composite);
        }
        let denominator = rem::<B>(&mul_u64::<B>(y1, 2)?, modulus)?;
        let x_squared = mul::<B>(x1, x1)?;
        let numerator = add::<B>(&mul_u64::<B>(&x_squared, 3)?, &curve.a)?;
        modular_mul::<B>(
            &numerator,
            &modular_inverse::<B>(&denominator, modulus)?,
            modulus,
        )?
    } else {
        let numerator = modular_sub::<B>(y2, y1, modulus)?;
        let denominator = modular_sub::<B>(x2, x1, modulus)?;
        modular_mul::<B>(
            &numerator,
            &modular_inverse::<B>(&denominator, modulus)?,
            modulus,
        )?
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
    for index in 0..crate::arithmetic::bit_length::<B>(scalar) {
        if bit::<B>(scalar, index) {
            output = point_add::<B>(curve, modulus, &output, &addend)?;
        }
        if index + 1 < crate::arithmetic::bit_length::<B>(scalar) {
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
}

impl Default for ProverOptions {
    fn default() -> Self {
        Self {
            trial_division_limit: 10_000,
            point_attempts: 128,
            max_depth: 64,
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
        let step = find_step::<B, R>(&current, rng, &primes, options.point_attempts)?;
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
}

fn find_step<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    candidate: &B,
    rng: &mut R,
    primes: &[u32],
    point_attempts: u32,
) -> Result<WorkingStep<B>> {
    let mut search = StepSearch::<B, R> {
        candidate,
        rng,
        primes,
        point_attempts,
    };
    let order = add_u64::<B>(candidate, 1)?;
    let modulo_four = rem::<B>(candidate, &from_u64::<B>(4)?)?;
    if cmp_u64::<B>(&modulo_four, 3)? == Ordering::Equal {
        let curve = AffineCurve::<B> {
            a: one::<B>()?,
            b: zero::<B>()?,
        };
        if let Some(step) = try_order::<B, R>(&mut search, &order, &curve, None)? {
            return Ok(step);
        }
    }
    let modulo_three = rem::<B>(candidate, &from_u64::<B>(3)?)?;
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
        let Some((trace, _)) = cornacchia::<B>(candidate, discriminant.value)? else {
            continue;
        };
        let Some(invariants) = j_invariants::<B>(candidate, discriminant.polynomial)? else {
            continue;
        };
        for invariant in invariants {
            let base = curve_from_j::<B>(candidate, &invariant)?;
            let twist = quadratic_twist::<B>(candidate, &base)?;
            let lower = sub::<B>(&order, &trace)?;
            let upper = add::<B>(&order, &trace)?;
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
    let Some(q) = split_order::<B, R>(search.candidate, order, search.primes, search.rng)? else {
        return Ok(None);
    };
    let cofactor = div::<B>(order, &q)?;
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
    rng: &mut R,
) -> Result<Option<B>> {
    let mut remaining = order.clone();
    for &prime in primes {
        let prime = from_u64::<B>(u64::from(prime))?;
        while is_zero::<B>(&rem::<B>(&remaining, &prime)?)? {
            remaining = div::<B>(&remaining, &prime)?;
        }
    }
    let fourth_root = ceil_fourth_root::<B>(candidate)?;
    let root_plus_one = add_u64::<B>(&fourth_root, 1)?;
    let bound = mul::<B>(&root_plus_one, &root_plus_one)?;
    find_large_prime_factor::<B, R>(&remaining, candidate, &bound, primes, rng, 0)
}

fn find_large_prime_factor<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    value: &B,
    candidate: &B,
    bound: &B,
    primes: &[u32],
    rng: &mut R,
    depth: usize,
) -> Result<Option<B>> {
    if value <= bound || depth > 32 {
        return Ok(None);
    }
    if is_probable_prime::<B>(value, primes)? {
        return Ok((value < candidate).then(|| value.clone()));
    }
    let Some(factor) = pollard_rho::<B, R>(value, rng)? else {
        return Ok(None);
    };
    let other = div::<B>(value, &factor)?;
    let (first, second) = if factor >= other {
        (factor, other)
    } else {
        (other, factor)
    };
    if let Some(found) =
        find_large_prime_factor::<B, R>(&first, candidate, bound, primes, rng, depth + 1)?
    {
        return Ok(Some(found));
    }
    find_large_prime_factor::<B, R>(&second, candidate, bound, primes, rng, depth + 1)
}

fn pollard_rho<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    value: &B,
    rng: &mut R,
) -> Result<Option<B>> {
    if is_even::<B>(value) {
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
        while divisor == one && iterations < MAX_ITERATIONS {
            x = y.clone();
            for _ in 0..power {
                y = rem::<B>(&add::<B>(&mul::<B>(&y, &y)?, &constant)?, value)?;
            }
            iterations += power;
            let mut offset = 0usize;
            while offset < power && divisor == one {
                saved_y = y.clone();
                let count = BATCH.min(power - offset);
                let mut product = one.clone();
                for _ in 0..count {
                    y = rem::<B>(&add::<B>(&mul::<B>(&y, &y)?, &constant)?, value)?;
                    let difference = if x >= y {
                        sub::<B>(&x, &y)?
                    } else {
                        sub::<B>(&y, &x)?
                    };
                    product = rem::<B>(&mul::<B>(&product, &difference)?, value)?;
                }
                divisor = gcd::<B>(&product, value)?;
                offset += count;
                iterations += count;
            }
            power = power.saturating_mul(2);
        }
        if divisor == *value {
            divisor = one.clone();
            while divisor == one && iterations < MAX_ITERATIONS * 2 {
                saved_y = rem::<B>(&add::<B>(&mul::<B>(&saved_y, &saved_y)?, &constant)?, value)?;
                let difference = if x >= saved_y {
                    sub::<B>(&x, &saved_y)?
                } else {
                    sub::<B>(&saved_y, &x)?
                };
                divisor = gcd::<B>(&difference, value)?;
                iterations += 1;
            }
        }
        if divisor != one && divisor != *value {
            return Ok(Some(divisor));
        }
    }
    Ok(None)
}

fn curve_from_j<B: ArithmeticBackend>(candidate: &B, invariant: &B) -> Result<AffineCurve<B>> {
    let denominator = modular_sub::<B>(&from_u64::<B>(1728)?, invariant, candidate)?;
    let inverse = modular_inverse::<B>(&denominator, candidate)?;
    let k = modular_mul::<B>(invariant, &inverse, candidate)?;
    Ok(AffineCurve {
        a: rem::<B>(&mul_u64::<B>(&k, 3)?, candidate)?,
        b: rem::<B>(&mul_u64::<B>(&k, 2)?, candidate)?,
    })
}

fn quadratic_twist<B: ArithmeticBackend>(
    candidate: &B,
    curve: &AffineCurve<B>,
) -> Result<AffineCurve<B>> {
    let mut non_residue = from_u64::<B>(2)?;
    while jacobi::<B>(&non_residue, candidate)? != -1 {
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
    let cofactor = div::<B>(order, q)?;
    for _ in 0..attempts {
        let x = random_below::<B, R>(modulus, rng)?;
        let x_squared = modular_mul::<B>(&x, &x, modulus)?;
        let x_cubed = modular_mul::<B>(&x_squared, &x, modulus)?;
        let ax = modular_mul::<B>(&curve.a, &x, modulus)?;
        let rhs = rem::<B>(&add::<B>(&add::<B>(&x_cubed, &ax)?, &curve.b)?, modulus)?;
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
    let byte_length = bit_length::<B>(modulus).div_ceil(8);
    let excess = byte_length * 8 - bit_length::<B>(modulus);
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
    if is_even::<B>(&n)
        || cmp_u64::<B>(&n, 3)? == Ordering::Less
        || q >= n
        || cmp_u64::<B>(&q, 2)? == Ordering::Less
    {
        return Err(Error::InvalidProof("invalid ECPP step integers"));
    }
    let fourth_root = ceil_fourth_root::<B>(&n)?;
    let root_plus_one = add_u64::<B>(&fourth_root, 1)?;
    if q <= mul::<B>(&root_plus_one, &root_plus_one)? {
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
    let discriminant = rem::<B>(
        &add::<B>(&mul_u64::<B>(&a_cubed, 4)?, &mul_u64::<B>(&b_squared, 27)?)?,
        &n,
    )?;
    if gcd::<B>(&discriminant, &n)? != one::<B>()? {
        return Err(Error::InvalidProof(
            "curve is singular modulo a divisor of n",
        ));
    }
    let x_squared = modular_mul::<B>(x, x, &n)?;
    let x_cubed = modular_mul::<B>(&x_squared, x, &n)?;
    let ax = modular_mul::<B>(&curve.a, x, &n)?;
    let rhs = rem::<B>(&add::<B>(&add::<B>(&x_cubed, &ax)?, &curve.b)?, &n)?;
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
