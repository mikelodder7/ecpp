//! Arithmetic backends used by the heap-backed ECPP engine.

use alloc::vec::Vec;
use core::ops::{Add, Div, Mul, Rem, Shl, Shr, Sub};

#[cfg(feature = "crypto-bigint")]
use crypto_bigint::{CheckedAdd, CheckedSub};

#[cfg(feature = "openssl")]
use core::cmp::Ordering;

use crate::{Error, Result};

/// A non-negative arbitrary-precision value usable by the generic ECPP engine.
///
/// Every arithmetic operator is fallible, so adapters such as OpenSSL can
/// propagate allocation or library failures. Implementations own the left-hand
/// operand and borrow the right-hand operand, which permits in-place arithmetic
/// without requiring the engine to copy both inputs.
///
/// Only the two byte codecs are required. The provided methods have correct
/// byte-based or operator-based defaults; backends should override them with
/// native implementations where the underlying library offers one.
pub trait ArithmeticBackend:
    Sized
    + Clone
    + Eq
    + Ord
    + for<'a> Add<&'a Self, Output = Result<Self>>
    + for<'a> Sub<&'a Self, Output = Result<Self>>
    + for<'a> Mul<&'a Self, Output = Result<Self>>
    + for<'a> Div<&'a Self, Output = Result<Self>>
    + for<'a> Rem<&'a Self, Output = Result<Self>>
    + Shl<usize, Output = Result<Self>>
    + Shr<usize, Output = Result<Self>>
{
    /// Decodes an unsigned big-endian magnitude.
    fn from_be_bytes(bytes: &[u8]) -> Result<Self>;

    /// Encodes this value without redundant leading zeros.
    fn to_be_bytes(&self) -> Vec<u8>;

    /// Returns the number of significant bits, so zero has bit length zero.
    fn bit_length(&self) -> usize {
        let bytes = self.to_be_bytes();
        bytes
            .first()
            .map_or(0, |first| bytes.len() * 8 - first.leading_zeros() as usize)
    }

    /// Returns the bit at `index`, where index zero is the least significant.
    ///
    /// Bits at or beyond the value's width are zero.
    fn bit(&self, index: usize) -> bool {
        let bytes = self.to_be_bytes();
        let byte_from_end = index / 8;
        if byte_from_end >= bytes.len() {
            return false;
        }
        bytes[bytes.len() - byte_from_end - 1] & (1 << (index % 8)) != 0
    }

    /// Reports whether this value is even.
    fn is_even(&self) -> bool {
        !self.bit(0)
    }

    /// Reports whether this value is zero.
    fn is_zero(&self) -> bool {
        self.bit_length() == 0
    }

    /// Reports whether this value is one.
    fn is_one(&self) -> bool {
        self.bit_length() == 1
    }

    /// Computes `self ^ exponent mod modulus`.
    fn modular_pow(&self, exponent: &Self, modulus: &Self) -> Result<Self> {
        generic_modular_pow(self, exponent, modulus)
    }

    /// Computes the greatest common divisor of `self` and `other`.
    fn gcd(&self, other: &Self) -> Result<Self> {
        let mut left = self.clone();
        let mut right = other.clone();
        while !right.is_zero() {
            let remainder = (left % &right)?;
            left = right;
            right = remainder;
        }
        Ok(left)
    }

    /// Computes the inverse of `self` modulo `modulus`.
    ///
    /// Returns [`Error::Composite`] when the inverse does not exist, matching
    /// the ECPP convention that a failed inversion witnesses compositeness.
    fn modular_inverse(&self, modulus: &Self) -> Result<Self> {
        let zero = Self::from_be_bytes(&[])?;
        let one = Self::from_be_bytes(&[1])?;
        let mut coefficient = zero.clone();
        let mut next_coefficient = one.clone();
        let mut remainder = modulus.clone();
        let mut next_remainder = (self.clone() % modulus)?;
        while !next_remainder.is_zero() {
            let quotient = (remainder.clone() / &next_remainder)?;
            let product = ((quotient.clone() * &next_coefficient)? % modulus)?;
            let old_coefficient = coefficient;
            coefficient = next_coefficient;
            next_coefficient = if old_coefficient >= product {
                ((old_coefficient - &product)? % modulus)?
            } else {
                let difference = ((product - &old_coefficient)? % modulus)?;
                if difference.is_zero() {
                    zero.clone()
                } else {
                    (modulus.clone() - &difference)?
                }
            };
            let following_remainder = (remainder - &(quotient * &next_remainder)?)?;
            remainder = next_remainder;
            next_remainder = following_remainder;
        }
        if !remainder.is_one() {
            return Err(Error::Composite);
        }
        Ok(coefficient)
    }

    /// Computes the Jacobi symbol `(self / modulus)`.
    ///
    /// Returns zero when `modulus` is zero or even.
    fn jacobi(&self, modulus: &Self) -> Result<i8> {
        if modulus.is_zero() || modulus.is_even() {
            return Ok(0);
        }
        let mut value = (self.clone() % modulus)?;
        let mut modulus = modulus.clone();
        let mut result = 1i8;
        while !value.is_zero() {
            while value.is_even() {
                value = (value >> 1)?;
                // For odd values, residue 3 or 5 modulo 8 is exactly when bits
                // one and two differ.
                if modulus.bit(1) != modulus.bit(2) {
                    result = -result;
                }
            }
            core::mem::swap(&mut value, &mut modulus);
            // Both operands are odd here, so bit one alone decides whether
            // each is 3 modulo 4.
            if value.bit(1) && modulus.bit(1) {
                result = -result;
            }
            value = (value % &modulus)?;
        }
        Ok(if modulus.is_one() { result } else { 0 })
    }
}

/// Square-and-multiply `base ^ exponent mod modulus` over the operator surface.
///
/// This is the [`ArithmeticBackend::modular_pow`] default; backends whose
/// native exponentiation has preconditions (such as an odd modulus) also use
/// it as their fallback.
pub(crate) fn generic_modular_pow<B: ArithmeticBackend>(
    base: &B,
    exponent: &B,
    modulus: &B,
) -> Result<B> {
    let mut output = B::from_be_bytes(&[1])?;
    let mut base = (base.clone() % modulus)?;
    let bits = exponent.bit_length();
    for index in 0..bits {
        if exponent.bit(index) {
            output = ((output * &base)? % modulus)?;
        }
        if index + 1 < bits {
            base = ((base.clone() * &base)? % modulus)?;
        }
    }
    Ok(output)
}

/// A `num-bigint` value for the backend-neutral engine.
#[cfg(feature = "num-bigint")]
#[cfg_attr(docsrs, doc(cfg(feature = "num-bigint")))]
#[derive(Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct NumBigint(num_bigint::BigUint);

#[cfg(feature = "num-bigint")]
impl ArithmeticBackend for NumBigint {
    fn from_be_bytes(bytes: &[u8]) -> Result<Self> {
        Ok(Self(num_bigint::BigUint::from_bytes_be(bytes)))
    }

    fn to_be_bytes(&self) -> Vec<u8> {
        self.0.to_bytes_be()
    }

    fn bit_length(&self) -> usize {
        self.0.bits() as usize
    }

    fn bit(&self, index: usize) -> bool {
        self.0.bit(index as u64)
    }

    fn modular_pow(&self, exponent: &Self, modulus: &Self) -> Result<Self> {
        // `BigUint::modpow` panics on a zero modulus.
        if modulus.0 == num_bigint::BigUint::default() {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        Ok(Self(self.0.modpow(&exponent.0, &modulus.0)))
    }

    fn modular_inverse(&self, modulus: &Self) -> Result<Self> {
        if modulus.0 == num_bigint::BigUint::default() {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        self.0.modinv(&modulus.0).map(Self).ok_or(Error::Composite)
    }
}

#[cfg(feature = "num-bigint")]
impl Add<&Self> for NumBigint {
    type Output = Result<Self>;

    fn add(mut self, right: &Self) -> Self::Output {
        self.0 += &right.0;
        Ok(self)
    }
}

#[cfg(feature = "num-bigint")]
impl Sub<&Self> for NumBigint {
    type Output = Result<Self>;

    fn sub(mut self, right: &Self) -> Self::Output {
        if self < *right {
            return Err(Error::Arithmetic("subtraction would be negative"));
        }
        self.0 -= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "num-bigint")]
impl Mul<&Self> for NumBigint {
    type Output = Result<Self>;

    fn mul(mut self, right: &Self) -> Self::Output {
        self.0 *= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "num-bigint")]
impl Div<&Self> for NumBigint {
    type Output = Result<Self>;

    fn div(mut self, right: &Self) -> Self::Output {
        if right.0 == num_bigint::BigUint::default() {
            return Err(Error::Arithmetic("division by zero"));
        }
        self.0 /= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "num-bigint")]
impl Rem<&Self> for NumBigint {
    type Output = Result<Self>;

    fn rem(mut self, right: &Self) -> Self::Output {
        if right.0 == num_bigint::BigUint::default() {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        self.0 %= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "num-bigint")]
impl Shl<usize> for NumBigint {
    type Output = Result<Self>;

    fn shl(mut self, bits: usize) -> Self::Output {
        self.0 <<= bits;
        Ok(self)
    }
}

#[cfg(feature = "num-bigint")]
impl Shr<usize> for NumBigint {
    type Output = Result<Self>;

    fn shr(mut self, bits: usize) -> Self::Output {
        self.0 >>= bits;
        Ok(self)
    }
}

/// A GMP-backed value through `rug`.
#[cfg(feature = "rug")]
#[cfg_attr(docsrs, doc(cfg(feature = "rug")))]
#[derive(Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct Rug(rug::Integer);

#[cfg(feature = "rug")]
impl ArithmeticBackend for Rug {
    fn from_be_bytes(bytes: &[u8]) -> Result<Self> {
        Ok(Self(rug::Integer::from_digits(
            bytes,
            rug::integer::Order::MsfBe,
        )))
    }

    fn to_be_bytes(&self) -> Vec<u8> {
        let mut bytes = alloc::vec![0u8; self.0.significant_digits::<u8>()];
        self.0.write_digits(&mut bytes, rug::integer::Order::MsfBe);
        bytes
    }

    fn bit_length(&self) -> usize {
        self.0.significant_bits() as usize
    }

    fn bit(&self, index: usize) -> bool {
        u32::try_from(index).is_ok_and(|index| self.0.get_bit(index))
    }

    fn modular_pow(&self, exponent: &Self, modulus: &Self) -> Result<Self> {
        // GMP exponentiation is undefined for a zero modulus.
        if modulus.0 == 0 {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        self.0
            .clone()
            .pow_mod(&exponent.0, &modulus.0)
            .map(Self)
            .map_err(|_| Error::Arithmetic("modular exponentiation failed"))
    }

    fn gcd(&self, other: &Self) -> Result<Self> {
        Ok(Self(self.0.clone().gcd(&other.0)))
    }

    fn modular_inverse(&self, modulus: &Self) -> Result<Self> {
        if modulus.0 == 0 {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        self.0
            .clone()
            .invert(&modulus.0)
            .map(Self)
            .map_err(|_| Error::Composite)
    }

    fn jacobi(&self, modulus: &Self) -> Result<i8> {
        // GMP defines the Jacobi symbol only for odd moduli.
        if modulus.0 == 0 || !modulus.0.get_bit(0) {
            return Ok(0);
        }
        Ok(rug::Integer::jacobi(&self.0, &modulus.0) as i8)
    }
}

#[cfg(feature = "rug")]
impl Add<&Self> for Rug {
    type Output = Result<Self>;

    fn add(mut self, right: &Self) -> Self::Output {
        self.0 += &right.0;
        Ok(self)
    }
}

#[cfg(feature = "rug")]
impl Sub<&Self> for Rug {
    type Output = Result<Self>;

    fn sub(mut self, right: &Self) -> Self::Output {
        if self < *right {
            return Err(Error::Arithmetic("subtraction would be negative"));
        }
        self.0 -= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "rug")]
impl Mul<&Self> for Rug {
    type Output = Result<Self>;

    fn mul(mut self, right: &Self) -> Self::Output {
        self.0 *= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "rug")]
impl Div<&Self> for Rug {
    type Output = Result<Self>;

    fn div(mut self, right: &Self) -> Self::Output {
        if right.0 == 0 {
            return Err(Error::Arithmetic("division by zero"));
        }
        self.0 /= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "rug")]
impl Rem<&Self> for Rug {
    type Output = Result<Self>;

    fn rem(mut self, right: &Self) -> Self::Output {
        if right.0 == 0 {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        self.0 %= &right.0;
        Ok(self)
    }
}

#[cfg(feature = "rug")]
impl Shl<usize> for Rug {
    type Output = Result<Self>;

    fn shl(mut self, bits: usize) -> Self::Output {
        self.0 <<= bits;
        Ok(self)
    }
}

#[cfg(feature = "rug")]
impl Shr<usize> for Rug {
    type Output = Result<Self>;

    fn shr(mut self, bits: usize) -> Self::Output {
        self.0 >>= bits;
        Ok(self)
    }
}

/// A heap-backed `crypto-bigint` value for the generic ECPP engine.
///
/// This adapter widens its operands before each operation, avoiding the
/// fixed-precision overflow behavior of a bare `BoxedUint`.
#[cfg(feature = "crypto-bigint")]
#[cfg_attr(docsrs, doc(cfg(all(feature = "crypto-bigint", feature = "alloc"))))]
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct CryptoBigint(crypto_bigint::BoxedUint);

#[cfg(feature = "crypto-bigint")]
impl CryptoBigint {
    fn normalize(value: crypto_bigint::BoxedUint) -> Self {
        let bytes = value.to_be_bytes_trimmed_vartime();
        if bytes.is_empty() {
            Self(crypto_bigint::BoxedUint::zero())
        } else {
            Self(crypto_bigint::BoxedUint::from_be_slice_vartime(&bytes))
        }
    }
}

#[cfg(feature = "crypto-bigint")]
impl Default for CryptoBigint {
    fn default() -> Self {
        Self(crypto_bigint::BoxedUint::zero())
    }
}

#[cfg(feature = "crypto-bigint")]
impl ArithmeticBackend for CryptoBigint {
    fn from_be_bytes(bytes: &[u8]) -> Result<Self> {
        if bytes.iter().all(|byte| *byte == 0) {
            Ok(Self::default())
        } else {
            Ok(Self(crypto_bigint::BoxedUint::from_be_slice_vartime(bytes)))
        }
    }

    fn to_be_bytes(&self) -> Vec<u8> {
        self.0.to_be_bytes_trimmed_vartime().into_vec()
    }

    fn bit_length(&self) -> usize {
        self.0.bits_vartime() as usize
    }

    fn bit(&self, index: usize) -> bool {
        u32::try_from(index).is_ok_and(|index| self.0.bit_vartime(index))
    }
}

#[cfg(feature = "crypto-bigint")]
impl Add<&Self> for CryptoBigint {
    type Output = Result<Self>;

    fn add(self, right: &Self) -> Self::Output {
        Ok(Self::normalize(self.0.concatenating_add(&right.0)))
    }
}

#[cfg(feature = "crypto-bigint")]
impl Sub<&Self> for CryptoBigint {
    type Output = Result<Self>;

    fn sub(self, right: &Self) -> Self::Output {
        if self < *right {
            return Err(Error::Arithmetic("subtraction would be negative"));
        }
        let (difference, _) = self.0.borrowing_sub(&right.0, crypto_bigint::Limb::ZERO);
        Ok(Self::normalize(difference))
    }
}

#[cfg(feature = "crypto-bigint")]
impl Mul<&Self> for CryptoBigint {
    type Output = Result<Self>;

    fn mul(self, right: &Self) -> Self::Output {
        let product = crypto_bigint::ConcatenatingMul::concatenating_mul(&self.0, &right.0);
        Ok(Self::normalize(product))
    }
}

#[cfg(feature = "crypto-bigint")]
impl Div<&Self> for CryptoBigint {
    type Output = Result<Self>;

    fn div(self, right: &Self) -> Self::Output {
        let Some(divisor) = Option::from(crypto_bigint::NonZero::new(right.0.clone())) else {
            return Err(Error::Arithmetic("division by zero"));
        };
        let (quotient, _) = self.0.div_rem_vartime(&divisor);
        Ok(Self::normalize(quotient))
    }
}

#[cfg(feature = "crypto-bigint")]
impl Rem<&Self> for CryptoBigint {
    type Output = Result<Self>;

    fn rem(self, right: &Self) -> Self::Output {
        let Some(divisor) = Option::from(crypto_bigint::NonZero::new(right.0.clone())) else {
            return Err(Error::Arithmetic("remainder by zero"));
        };
        let (_, remainder) = self.0.div_rem_vartime(&divisor);
        Ok(Self::normalize(remainder))
    }
}

#[cfg(feature = "crypto-bigint")]
impl Shl<usize> for CryptoBigint {
    type Output = Result<Self>;

    fn shl(self, bits: usize) -> Self::Output {
        use crypto_bigint::Resize;

        let bits = u32::try_from(bits).map_err(|_| Error::Arithmetic("shift is too large"))?;
        if self.0.bits() == 0 {
            return Ok(Self::default());
        }
        let precision = self
            .0
            .bits()
            .checked_add(bits)
            .ok_or(Error::Arithmetic("shift is too large"))?
            .max(1)
            .next_multiple_of(crypto_bigint::Limb::BITS);
        let mut value = self.0.resize_unchecked(precision);
        value.shl_assign(bits);
        Ok(Self::normalize(value))
    }
}

#[cfg(feature = "crypto-bigint")]
impl Shr<usize> for CryptoBigint {
    type Output = Result<Self>;

    fn shr(self, bits: usize) -> Self::Output {
        let bits = u32::try_from(bits).map_err(|_| Error::Arithmetic("shift is too large"))?;
        if bits >= self.0.bits_precision() {
            return Ok(Self::default());
        }
        let mut value = self.0;
        value.shr_assign(bits);
        Ok(Self::normalize(value))
    }
}

/// A fixed-width `crypto-bigint::Uint` value for the generic ECPP engine.
///
/// `LIMBS` is the arithmetic working precision, not necessarily the width of
/// the candidate supplied to `prime::prove_with_backend`. It should normally
/// be at least twice the candidate width so unreduced products fit. Operations
/// return `Error::Arithmetic` when that working precision is exhausted.
#[cfg(feature = "crypto-bigint")]
#[cfg_attr(docsrs, doc(cfg(all(feature = "crypto-bigint", feature = "alloc"))))]
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct CryptoUint<const LIMBS: usize>(crypto_bigint::Uint<LIMBS>);

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> ArithmeticBackend for CryptoUint<LIMBS> {
    fn from_be_bytes(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > crypto_bigint::Uint::<LIMBS>::BYTES {
            return Err(Error::Arithmetic(
                "integer exceeds crypto-bigint working precision",
            ));
        }
        let mut padded = alloc::vec![0u8; crypto_bigint::Uint::<LIMBS>::BYTES];
        let start = padded.len() - bytes.len();
        padded[start..].copy_from_slice(bytes);
        Ok(Self(crypto_bigint::Uint::from_be_slice(&padded)))
    }

    fn to_be_bytes(&self) -> Vec<u8> {
        let bytes = self.0.to_be_bytes();
        let first = bytes
            .iter()
            .position(|byte| *byte != 0)
            .unwrap_or(bytes.len());
        bytes[first..].to_vec()
    }

    fn bit_length(&self) -> usize {
        self.0.bits_vartime() as usize
    }

    fn bit(&self, index: usize) -> bool {
        u32::try_from(index).is_ok_and(|index| {
            index < crypto_bigint::Uint::<LIMBS>::BITS && self.0.bit_vartime(index)
        })
    }

    fn modular_pow(&self, exponent: &Self, modulus: &Self) -> Result<Self> {
        use crypto_bigint::modular::{FixedMontyForm, FixedMontyParams};

        let Some(odd) = crypto_bigint::Odd::new(modulus.0).into_option() else {
            // Montgomery form requires an odd modulus.
            return generic_modular_pow(self, exponent, modulus);
        };
        let params = FixedMontyParams::new_vartime(odd);
        let reduced = self.0.rem_vartime(odd.as_nz_ref());
        let value = FixedMontyForm::new(&reduced, &params);
        Ok(Self(value.pow_vartime(&exponent.0).retrieve()))
    }

    fn gcd(&self, other: &Self) -> Result<Self> {
        Ok(Self(self.0.gcd_vartime(&other.0)))
    }

    fn modular_inverse(&self, modulus: &Self) -> Result<Self> {
        let Some(modulus) = crypto_bigint::NonZero::new(modulus.0).into_option() else {
            return Err(Error::Arithmetic("remainder by zero"));
        };
        self.0
            .rem_vartime(&modulus)
            .invert_mod(&modulus)
            .into_option()
            .map(Self)
            .ok_or(Error::Composite)
    }

    fn jacobi(&self, modulus: &Self) -> Result<i8> {
        let Some(odd) = crypto_bigint::Odd::new(modulus.0).into_option() else {
            return Ok(0);
        };
        let reduced = self.0.rem_vartime(odd.as_nz_ref());
        Ok(reduced.jacobi_symbol_vartime(&odd).into())
    }
}

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> Add<&Self> for CryptoUint<LIMBS> {
    type Output = Result<Self>;

    fn add(self, right: &Self) -> Self::Output {
        self.0
            .checked_add(&right.0)
            .into_option()
            .map(Self)
            .ok_or(Error::Arithmetic(
                "crypto-bigint addition exceeded working precision",
            ))
    }
}

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> Sub<&Self> for CryptoUint<LIMBS> {
    type Output = Result<Self>;

    fn sub(self, right: &Self) -> Self::Output {
        self.0
            .checked_sub(&right.0)
            .into_option()
            .map(Self)
            .ok_or(Error::Arithmetic("subtraction would be negative"))
    }
}

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> Mul<&Self> for CryptoUint<LIMBS> {
    type Output = Result<Self>;

    fn mul(self, right: &Self) -> Self::Output {
        self.0
            .checked_mul(&right.0)
            .into_option()
            .map(Self)
            .ok_or(Error::Arithmetic(
                "crypto-bigint multiplication exceeded working precision",
            ))
    }
}

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> Div<&Self> for CryptoUint<LIMBS> {
    type Output = Result<Self>;

    fn div(self, right: &Self) -> Self::Output {
        let Some(divisor) = crypto_bigint::NonZero::new(right.0).into_option() else {
            return Err(Error::Arithmetic("division by zero"));
        };
        let (quotient, _) = self.0.div_rem_vartime(&divisor);
        Ok(Self(quotient))
    }
}

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> Rem<&Self> for CryptoUint<LIMBS> {
    type Output = Result<Self>;

    fn rem(self, right: &Self) -> Self::Output {
        let Some(divisor) = crypto_bigint::NonZero::new(right.0).into_option() else {
            return Err(Error::Arithmetic("remainder by zero"));
        };
        let (_, remainder) = self.0.div_rem_vartime(&divisor);
        Ok(Self(remainder))
    }
}

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> Shl<usize> for CryptoUint<LIMBS> {
    type Output = Result<Self>;

    fn shl(self, bits: usize) -> Self::Output {
        let bits = u32::try_from(bits).map_err(|_| Error::Arithmetic("shift is too large"))?;
        if self
            .0
            .bits_vartime()
            .checked_add(bits)
            .is_none_or(|width| width > crypto_bigint::Uint::<LIMBS>::BITS)
        {
            return Err(Error::Arithmetic(
                "crypto-bigint shift exceeded working precision",
            ));
        }
        Ok(Self(self.0 << bits))
    }
}

#[cfg(feature = "crypto-bigint")]
impl<const LIMBS: usize> Shr<usize> for CryptoUint<LIMBS> {
    type Output = Result<Self>;

    fn shr(self, bits: usize) -> Self::Output {
        let bits = u32::try_from(bits).map_err(|_| Error::Arithmetic("shift is too large"))?;
        if bits >= crypto_bigint::Uint::<LIMBS>::BITS {
            Ok(Self(crypto_bigint::Uint::ZERO))
        } else {
            Ok(Self(self.0 >> bits))
        }
    }
}

/// An OpenSSL `BIGNUM` value.
///
/// The value is held as canonical bytes between operations. This makes cloning
/// infallible and keeps OpenSSL allocation failures in ordinary `Result` paths.
#[cfg(feature = "openssl")]
#[cfg_attr(docsrs, doc(cfg(feature = "openssl")))]
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct OpenSsl(Vec<u8>);

#[cfg(feature = "openssl")]
impl Ord for OpenSsl {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0
            .len()
            .cmp(&other.0.len())
            .then_with(|| self.0.cmp(&other.0))
    }
}

#[cfg(feature = "openssl")]
impl PartialOrd for OpenSsl {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(feature = "openssl")]
impl OpenSsl {
    fn decode(bytes: &[u8]) -> Result<openssl::bn::BigNum> {
        let encoded = if bytes.is_empty() { &[0u8][..] } else { bytes };
        openssl::bn::BigNum::from_slice(encoded)
            .map_err(|_| Error::Arithmetic("OpenSSL could not allocate a BIGNUM"))
    }

    fn encode(value: &openssl::bn::BigNumRef) -> Vec<u8> {
        canonicalize(value.to_vec())
    }

    fn binary(
        &self,
        right: &Self,
        operation: impl FnOnce(
            &mut openssl::bn::BigNum,
            &openssl::bn::BigNumRef,
            &openssl::bn::BigNumRef,
            &mut openssl::bn::BigNumContextRef,
        ) -> core::result::Result<(), openssl::error::ErrorStack>,
    ) -> Result<Self> {
        let left = Self::decode(&self.0)?;
        let right = Self::decode(&right.0)?;
        let mut output = openssl::bn::BigNum::new()
            .map_err(|_| Error::Arithmetic("OpenSSL could not allocate a BIGNUM"))?;
        let mut context = openssl::bn::BigNumContext::new()
            .map_err(|_| Error::Arithmetic("OpenSSL could not allocate a BIGNUM context"))?;
        operation(&mut output, &left, &right, &mut context)
            .map_err(|_| Error::Arithmetic("OpenSSL BIGNUM arithmetic failed"))?;
        Ok(Self(Self::encode(&output)))
    }

    fn shift(self, bits: usize, left: bool) -> Result<Self> {
        let bits = i32::try_from(bits).map_err(|_| Error::Arithmetic("shift is too large"))?;
        let value = Self::decode(&self.0)?;
        let mut output = openssl::bn::BigNum::new()
            .map_err(|_| Error::Arithmetic("OpenSSL could not allocate a BIGNUM"))?;
        if left {
            output
                .lshift(&value, bits)
                .map_err(|_| Error::Arithmetic("OpenSSL BIGNUM shift failed"))?;
        } else {
            output
                .rshift(&value, bits)
                .map_err(|_| Error::Arithmetic("OpenSSL BIGNUM shift failed"))?;
        }
        Ok(Self(Self::encode(&output)))
    }
}

#[cfg(feature = "openssl")]
impl ArithmeticBackend for OpenSsl {
    fn from_be_bytes(bytes: &[u8]) -> Result<Self> {
        Ok(Self(canonicalize(bytes.to_vec())))
    }

    fn to_be_bytes(&self) -> Vec<u8> {
        self.0.clone()
    }

    fn bit_length(&self) -> usize {
        self.0
            .first()
            .map_or(0, |first| self.0.len() * 8 - first.leading_zeros() as usize)
    }

    fn bit(&self, index: usize) -> bool {
        let byte_from_end = index / 8;
        if byte_from_end >= self.0.len() {
            return false;
        }
        self.0[self.0.len() - byte_from_end - 1] & (1 << (index % 8)) != 0
    }

    fn modular_pow(&self, exponent: &Self, modulus: &Self) -> Result<Self> {
        if modulus.0.is_empty() {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        let base = Self::decode(&self.0)?;
        let exponent = Self::decode(&exponent.0)?;
        let modulus = Self::decode(&modulus.0)?;
        let mut output = openssl::bn::BigNum::new()
            .map_err(|_| Error::Arithmetic("OpenSSL could not allocate a BIGNUM"))?;
        let mut context = openssl::bn::BigNumContext::new()
            .map_err(|_| Error::Arithmetic("OpenSSL could not allocate a BIGNUM context"))?;
        output
            .mod_exp(&base, &exponent, &modulus, &mut context)
            .map_err(|_| Error::Arithmetic("OpenSSL BIGNUM arithmetic failed"))?;
        Ok(Self(Self::encode(&output)))
    }

    fn gcd(&self, other: &Self) -> Result<Self> {
        self.binary(other, |output, left, right, context| {
            openssl::bn::BigNumRef::gcd(output, left, right, context)
        })
    }

    fn modular_inverse(&self, modulus: &Self) -> Result<Self> {
        if modulus.0.is_empty() {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        // `BN_mod_inverse` reports "no inverse" and real failures on the same
        // error path, so witness compositeness with an explicit gcd first.
        if self.gcd(modulus)? != Self::from_be_bytes(&[1])? {
            return Err(Error::Composite);
        }
        self.binary(modulus, |output, value, modulus, context| {
            openssl::bn::BigNumRef::mod_inverse(output, value, modulus, context)
        })
    }
}

#[cfg(feature = "openssl")]
impl Add<&Self> for OpenSsl {
    type Output = Result<Self>;

    fn add(self, right: &Self) -> Self::Output {
        self.binary(right, |output, left, right, _| {
            output.checked_add(left, right)
        })
    }
}

#[cfg(feature = "openssl")]
impl Sub<&Self> for OpenSsl {
    type Output = Result<Self>;

    fn sub(self, right: &Self) -> Self::Output {
        if self < *right {
            return Err(Error::Arithmetic("subtraction would be negative"));
        }
        self.binary(right, |output, left, right, _| {
            output.checked_sub(left, right)
        })
    }
}

#[cfg(feature = "openssl")]
impl Mul<&Self> for OpenSsl {
    type Output = Result<Self>;

    fn mul(self, right: &Self) -> Self::Output {
        self.binary(right, |output, left, right, context| {
            output.checked_mul(left, right, context)
        })
    }
}

#[cfg(feature = "openssl")]
impl Div<&Self> for OpenSsl {
    type Output = Result<Self>;

    fn div(self, right: &Self) -> Self::Output {
        if right.0.is_empty() {
            return Err(Error::Arithmetic("division by zero"));
        }
        self.binary(right, |output, left, right, context| {
            output.checked_div(left, right, context)
        })
    }
}

#[cfg(feature = "openssl")]
impl Rem<&Self> for OpenSsl {
    type Output = Result<Self>;

    fn rem(self, right: &Self) -> Self::Output {
        if right.0.is_empty() {
            return Err(Error::Arithmetic("remainder by zero"));
        }
        self.binary(right, |output, left, right, context| {
            output.checked_rem(left, right, context)
        })
    }
}

#[cfg(feature = "openssl")]
impl Shl<usize> for OpenSsl {
    type Output = Result<Self>;

    fn shl(self, bits: usize) -> Self::Output {
        self.shift(bits, true)
    }
}

#[cfg(feature = "openssl")]
impl Shr<usize> for OpenSsl {
    type Output = Result<Self>;

    fn shr(self, bits: usize) -> Self::Output {
        self.shift(bits, false)
    }
}

#[cfg(feature = "openssl")]
fn canonicalize(bytes: Vec<u8>) -> Vec<u8> {
    let first = bytes
        .iter()
        .position(|byte| *byte != 0)
        .unwrap_or(bytes.len());
    bytes[first..].to_vec()
}

pub(crate) fn from_u64<B: ArithmeticBackend>(value: u64) -> Result<B> {
    B::from_be_bytes(&value.to_be_bytes())
}

pub(crate) fn from_u128<B: ArithmeticBackend>(value: u128) -> Result<B> {
    B::from_be_bytes(&value.to_be_bytes())
}

pub(crate) fn to_u64<B: ArithmeticBackend>(value: &B) -> Option<u64> {
    let bytes = value.to_be_bytes();
    if bytes.len() > u64::BITS as usize / 8 {
        return None;
    }
    let mut output = [0u8; u64::BITS as usize / 8];
    let start = output.len() - bytes.len();
    output[start..].copy_from_slice(&bytes);
    Some(u64::from_be_bytes(output))
}

pub(crate) fn cmp_u64<B: ArithmeticBackend>(value: &B, other: u64) -> Result<core::cmp::Ordering> {
    Ok(value.cmp(&from_u64::<B>(other)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn operator_round_trip<B: ArithmeticBackend>() {
        let seven = B::from_be_bytes(&[7]).unwrap();
        let three = B::from_be_bytes(&[3]).unwrap();
        assert_eq!((seven.clone() + &three).unwrap().to_be_bytes(), [10]);
        assert_eq!((seven.clone() - &three).unwrap().to_be_bytes(), [4]);
        assert_eq!((seven.clone() * &three).unwrap().to_be_bytes(), [21]);
        assert_eq!((seven.clone() / &three).unwrap().to_be_bytes(), [2]);
        assert_eq!((seven.clone() % &three).unwrap().to_be_bytes(), [1]);
        assert_eq!((seven.clone() << 2).unwrap().to_be_bytes(), [28]);
        assert_eq!((seven >> 1).unwrap().to_be_bytes(), [3]);

        let zero = B::from_be_bytes(&[]).unwrap();
        assert!((three.clone() - &B::from_be_bytes(&[7]).unwrap()).is_err());
        assert!((three.clone() / &zero).is_err());
        assert!((three % &zero).is_err());
    }

    fn method_round_trip<B: ArithmeticBackend>() {
        let zero = B::from_be_bytes(&[]).unwrap();
        let one = B::from_be_bytes(&[1]).unwrap();
        let two = B::from_be_bytes(&[2]).unwrap();
        let three = B::from_be_bytes(&[3]).unwrap();
        let four = B::from_be_bytes(&[4]).unwrap();
        let five = B::from_be_bytes(&[5]).unwrap();
        let seven = B::from_be_bytes(&[7]).unwrap();
        let eight = B::from_be_bytes(&[8]).unwrap();

        assert_eq!(zero.bit_length(), 0);
        assert_eq!(one.bit_length(), 1);
        assert_eq!(seven.bit_length(), 3);
        assert_eq!(B::from_be_bytes(&[1, 0]).unwrap().bit_length(), 9);
        assert!(seven.bit(0) && seven.bit(1) && seven.bit(2));
        assert!(!seven.bit(3));
        assert!(!seven.bit(1_000));
        assert!(!seven.is_even());
        assert!(eight.is_even());
        assert!(zero.is_even());
        assert!(zero.is_zero());
        assert!(!one.is_zero());
        assert!(one.is_one());
        assert!(!two.is_one());
        assert!(!zero.is_one());

        assert_eq!(seven.modular_pow(&three, &five).unwrap().to_be_bytes(), [3]);
        assert_eq!(seven.modular_pow(&zero, &five).unwrap().to_be_bytes(), [1]);
        assert_eq!(three.modular_pow(&four, &eight).unwrap().to_be_bytes(), [1]);
        assert!(seven.modular_pow(&three, &zero).is_err());

        let twelve = B::from_be_bytes(&[12]).unwrap();
        let eighteen = B::from_be_bytes(&[18]).unwrap();
        assert_eq!(twelve.gcd(&eighteen).unwrap().to_be_bytes(), [6]);
        assert_eq!(zero.gcd(&seven).unwrap().to_be_bytes(), [7]);
        assert_eq!(seven.gcd(&zero).unwrap().to_be_bytes(), [7]);

        assert_eq!(three.modular_inverse(&seven).unwrap().to_be_bytes(), [5]);
        assert!(matches!(
            two.modular_inverse(&four),
            Err(crate::Error::Composite)
        ));
        assert!(matches!(
            zero.modular_inverse(&seven),
            Err(crate::Error::Composite)
        ));
        assert!(two.modular_inverse(&zero).is_err());

        assert_eq!(two.jacobi(&seven).unwrap(), 1);
        assert_eq!(three.jacobi(&seven).unwrap(), -1);
        assert_eq!(two.jacobi(&B::from_be_bytes(&[15]).unwrap()).unwrap(), 1);
        assert_eq!(seven.jacobi(&eight).unwrap(), 0);
        assert_eq!(seven.jacobi(&zero).unwrap(), 0);
        assert_eq!(zero.jacobi(&seven).unwrap(), 0);
    }

    #[test]
    #[cfg(feature = "num-bigint")]
    fn num_bigint_operators() {
        operator_round_trip::<NumBigint>();
        method_round_trip::<NumBigint>();
    }

    #[test]
    #[cfg(feature = "rug")]
    fn rug_operators() {
        operator_round_trip::<Rug>();
        method_round_trip::<Rug>();
    }

    #[test]
    #[cfg(feature = "crypto-bigint")]
    fn crypto_bigint_operators() {
        operator_round_trip::<CryptoBigint>();
        operator_round_trip::<CryptoUint<8>>();
        method_round_trip::<CryptoBigint>();
        method_round_trip::<CryptoUint<8>>();
    }

    #[test]
    #[cfg(feature = "openssl")]
    fn openssl_operators() {
        operator_round_trip::<OpenSsl>();
        method_round_trip::<OpenSsl>();
    }
}
