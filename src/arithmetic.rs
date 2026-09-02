//! Arithmetic backends used by the heap-backed ECPP engine.

use alloc::vec::Vec;
use core::ops::{Add, Div, Mul, Rem, Shl, Shr, Sub};

#[cfg(feature = "crypto-bigint")]
use crypto_bigint::{CheckedAdd, CheckedSub};

#[cfg(feature = "openssl")]
use core::cmp::Ordering;

#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
use crate::Error;
use crate::Result;

/// A non-negative arbitrary-precision value usable by the generic ECPP engine.
///
/// Every arithmetic operator is fallible, so adapters such as OpenSSL can
/// propagate allocation or library failures. Implementations own the left-hand
/// operand and borrow the right-hand operand, which permits in-place arithmetic
/// without requiring the engine to copy both inputs.
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
        self,
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

pub(crate) fn bit_length<B: ArithmeticBackend>(value: &B) -> usize {
    let bytes = value.to_be_bytes();
    bytes
        .first()
        .map_or(0, |first| bytes.len() * 8 - first.leading_zeros() as usize)
}

pub(crate) fn bit<B: ArithmeticBackend>(value: &B, index: usize) -> bool {
    let bytes = value.to_be_bytes();
    let byte_from_end = index / 8;
    if byte_from_end >= bytes.len() {
        return false;
    }
    bytes[bytes.len() - byte_from_end - 1] & (1 << (index % 8)) != 0
}

pub(crate) fn is_even<B: ArithmeticBackend>(value: &B) -> bool {
    !bit::<B>(value, 0)
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

    #[test]
    #[cfg(feature = "num-bigint")]
    fn num_bigint_operators() {
        operator_round_trip::<NumBigint>();
    }

    #[test]
    #[cfg(feature = "rug")]
    fn rug_operators() {
        operator_round_trip::<Rug>();
    }

    #[test]
    #[cfg(feature = "crypto-bigint")]
    fn crypto_bigint_operators() {
        operator_round_trip::<CryptoBigint>();
        operator_round_trip::<CryptoUint<8>>();
    }

    #[test]
    #[cfg(feature = "openssl")]
    fn openssl_operators() {
        operator_round_trip::<OpenSsl>();
    }
}
