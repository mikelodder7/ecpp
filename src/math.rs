use alloc::vec::Vec;

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer as _;
use num_traits::{One, ToPrimitive, Zero};

use crate::{Error, Result};

pub(crate) const SMALL_PRIME_LIMIT: u32 = 10_000;

pub(crate) fn mod_sub(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        (left - right) % modulus
    } else {
        (modulus - ((right - left) % modulus)) % modulus
    }
}

pub(crate) fn mod_signed(value: i128, modulus: &BigUint) -> BigUint {
    if value >= 0 {
        BigUint::from(value as u128) % modulus
    } else {
        let magnitude = BigUint::from(value.unsigned_abs()) % modulus;
        if magnitude.is_zero() {
            magnitude
        } else {
            modulus - magnitude
        }
    }
}

pub(crate) fn mod_inverse(value: &BigUint, modulus: &BigUint) -> Result<BigUint> {
    let value = value % modulus;
    let value_i = BigInt::from_biguint(Sign::Plus, value);
    let modulus_i = BigInt::from_biguint(Sign::Plus, modulus.clone());
    let extended = value_i.extended_gcd(&modulus_i);
    if extended.gcd != BigInt::one() {
        return Err(Error::Composite);
    }
    let mut inverse = extended.x % &modulus_i;
    if inverse.sign() == Sign::Minus {
        inverse += &modulus_i;
    }
    inverse.to_biguint().ok_or(Error::Composite)
}

pub(crate) fn integer_sqrt(value: &BigUint) -> BigUint {
    value.sqrt()
}

pub(crate) fn is_square(value: &BigUint) -> Option<BigUint> {
    let root = integer_sqrt(value);
    if &root * &root == *value {
        Some(root)
    } else {
        None
    }
}

pub(crate) fn jacobi(a: &BigUint, n: &BigUint) -> i8 {
    if n.is_zero() || n.is_even() {
        return 0;
    }
    let mut a = a % n;
    let mut n = n.clone();
    let mut result = 1i8;
    while !a.is_zero() {
        while a.is_even() {
            a >>= 1usize;
            let residue = (&n % 8u8).to_u8().unwrap_or(0);
            if residue == 3 || residue == 5 {
                result = -result;
            }
        }
        core::mem::swap(&mut a, &mut n);
        if (&a % 4u8 == BigUint::from(3u8)) && (&n % 4u8 == BigUint::from(3u8)) {
            result = -result;
        }
        a %= &n;
    }
    if n.is_one() { result } else { 0 }
}

/// Tonelli–Shanks. The caller only uses this after probable-prime screening.
pub(crate) fn modular_sqrt(value: &BigUint, modulus: &BigUint) -> Option<BigUint> {
    let value = value % modulus;
    if value.is_zero() {
        return Some(BigUint::zero());
    }
    if modulus == &BigUint::from(2u8) {
        return Some(value);
    }
    if jacobi(&value, modulus) != 1 {
        return None;
    }
    if modulus % 4u8 == BigUint::from(3u8) {
        return Some(value.modpow(&((modulus + 1u8) >> 2usize), modulus));
    }

    let mut odd = modulus - 1u8;
    let mut exponent = 0u32;
    while odd.is_even() {
        odd >>= 1usize;
        exponent += 1;
    }

    let mut non_residue = BigUint::from(2u8);
    while jacobi(&non_residue, modulus) != -1 {
        non_residue += 1u8;
        if &non_residue >= modulus {
            return None;
        }
    }

    let mut c = non_residue.modpow(&odd, modulus);
    let mut x = value.modpow(&((&odd + 1u8) >> 1usize), modulus);
    let mut t = value.modpow(&odd, modulus);
    let mut m = exponent;
    while !t.is_one() {
        let mut i = 1u32;
        let mut power = (&t * &t) % modulus;
        while !power.is_one() {
            power = (&power * &power) % modulus;
            i += 1;
            if i >= m {
                return None;
            }
        }
        let shift = (m - i - 1) as usize;
        let b = c.modpow(&(BigUint::one() << shift), modulus);
        x = (&x * &b) % modulus;
        let b_squared = (&b * &b) % modulus;
        t = (&t * &b_squared) % modulus;
        c = b_squared;
        m = i;
    }
    Some(x)
}

pub(crate) fn small_primes(limit: u32) -> Vec<u32> {
    let mut composite = alloc::vec![false; limit as usize + 1];
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

pub(crate) fn is_probable_prime(candidate: &BigUint, primes: &[u32]) -> bool {
    if candidate < &BigUint::from(2u8) {
        return false;
    }
    for &prime in primes {
        let prime_big = BigUint::from(prime);
        if candidate == &prime_big {
            return true;
        }
        if candidate % prime == BigUint::zero() {
            return false;
        }
    }
    if candidate.is_even() {
        return false;
    }

    const BASES: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    BASES
        .iter()
        .all(|base| miller_rabin(candidate, &BigUint::from(*base)))
}

fn miller_rabin(candidate: &BigUint, base: &BigUint) -> bool {
    if base >= candidate {
        return true;
    }
    let one = BigUint::one();
    let minus_one = candidate - &one;
    let mut odd = minus_one.clone();
    let mut exponent = 0u32;
    while odd.is_even() {
        odd >>= 1usize;
        exponent += 1;
    }
    let mut value = base.modpow(&odd, candidate);
    if value == one || value == minus_one {
        return true;
    }
    for _ in 1..exponent {
        value = (&value * &value) % candidate;
        if value == minus_one {
            return true;
        }
        if value == one {
            return false;
        }
    }
    false
}

pub(crate) fn is_prime_u64(candidate: u64) -> bool {
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
    // Deterministic for every 64-bit input.
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
