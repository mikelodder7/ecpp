use alloc::vec::Vec;
use core::fmt;

#[cfg(any(feature = "crypto-bigint", feature = "rug"))]
use alloc::vec;

use num_bigint::{BigInt, BigUint, Sign};

/// Conversion boundary between ECPP and an integer backend.
///
/// Primality is public-data computation, so canonical variable-time encoding is
/// intentional. Implementing this trait is sufficient to use another backend.
pub trait Integer: Sized {
    /// Returns the unsigned, big-endian magnitude with no redundant leading zero.
    fn to_be_bytes(&self) -> Vec<u8>;

    /// Constructs an integer from an unsigned, big-endian magnitude.
    ///
    /// Returns `None` if the magnitude does not fit the destination type.
    fn from_be_bytes(bytes: &[u8]) -> Option<Self>;

    /// Reports whether this integer is negative.
    fn is_negative(&self) -> bool {
        false
    }
}

/// A canonical, backend-neutral non-negative integer used in certificates.
#[derive(Clone, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Natural(Vec<u8>);

impl Natural {
    /// Creates a canonical value from big-endian bytes.
    pub fn from_be_bytes(bytes: &[u8]) -> Self {
        let first = bytes
            .iter()
            .position(|byte| *byte != 0)
            .unwrap_or(bytes.len());
        Self(bytes[first..].to_vec())
    }

    /// Returns the canonical big-endian magnitude.
    pub fn as_be_bytes(&self) -> &[u8] {
        &self.0
    }

    /// Converts this value into a supported integer backend.
    pub fn to_integer<T: Integer>(&self) -> Option<T> {
        T::from_be_bytes(&self.0)
    }

    pub(crate) fn from_biguint(value: &BigUint) -> Self {
        Self::from_be_bytes(&value.to_bytes_be())
    }

    pub(crate) fn to_biguint(&self) -> BigUint {
        BigUint::from_bytes_be(&self.0)
    }
}

impl fmt::Debug for Natural {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "0x")?;
        if self.0.is_empty() {
            return write!(formatter, "0");
        }
        for byte in &self.0 {
            write!(formatter, "{byte:02x}")?;
        }
        Ok(())
    }
}

impl fmt::Display for Natural {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}", self.to_biguint())
    }
}

impl Integer for BigUint {
    fn to_be_bytes(&self) -> Vec<u8> {
        self.to_bytes_be()
    }

    fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
        Some(Self::from_bytes_be(bytes))
    }
}

impl Integer for BigInt {
    fn to_be_bytes(&self) -> Vec<u8> {
        self.magnitude().to_bytes_be()
    }

    fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
        Some(Self::from_bytes_be(Sign::Plus, bytes))
    }

    fn is_negative(&self) -> bool {
        self.sign() == Sign::Minus
    }
}

impl Integer for Natural {
    fn to_be_bytes(&self) -> Vec<u8> {
        self.0.clone()
    }

    fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
        Some(Self::from_be_bytes(bytes))
    }
}

macro_rules! primitive_integer {
    ($($type:ty),+ $(,)?) => {$(
        impl Integer for $type {
            fn to_be_bytes(&self) -> Vec<u8> {
                let bytes = <$type>::to_be_bytes(*self);
                let first = bytes.iter().position(|byte| *byte != 0).unwrap_or(bytes.len());
                bytes[first..].to_vec()
            }

            fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
                if bytes.len() > core::mem::size_of::<$type>() {
                    return None;
                }
                let mut output = [0u8; core::mem::size_of::<$type>()];
                output[core::mem::size_of::<$type>() - bytes.len()..].copy_from_slice(bytes);
                Some(<$type>::from_be_bytes(output))
            }
        }
    )+};
}

primitive_integer!(u64, u128);

#[cfg(feature = "crypto-bigint")]
#[cfg_attr(docsrs, doc(cfg(feature = "crypto-bigint")))]
impl<const LIMBS: usize> Integer for crypto_bigint::Uint<LIMBS> {
    fn to_be_bytes(&self) -> Vec<u8> {
        let bytes = crypto_bigint::Uint::<LIMBS>::to_be_bytes(self);
        let first = bytes
            .iter()
            .position(|byte| *byte != 0)
            .unwrap_or(bytes.len());
        bytes[first..].to_vec()
    }

    fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
        let length = crypto_bigint::Limb::BYTES * LIMBS;
        if bytes.len() > length {
            return None;
        }
        let mut padded = vec![0u8; length];
        padded[length - bytes.len()..].copy_from_slice(bytes);
        Some(Self::from_be_slice(&padded))
    }
}

#[cfg(feature = "crypto-bigint")]
#[cfg_attr(docsrs, doc(cfg(feature = "crypto-bigint")))]
impl Integer for crypto_bigint::BoxedUint {
    fn to_be_bytes(&self) -> Vec<u8> {
        self.to_be_bytes_trimmed_vartime().into_vec()
    }

    fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
        Some(Self::from_be_slice_vartime(bytes))
    }
}

#[cfg(feature = "rug")]
#[cfg_attr(docsrs, doc(cfg(feature = "rug")))]
impl Integer for rug::Integer {
    fn to_be_bytes(&self) -> Vec<u8> {
        let mut bytes = vec![0u8; self.significant_digits::<u8>()];
        self.write_digits(&mut bytes, rug::integer::Order::MsfBe);
        bytes
    }

    fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
        Some(Self::from_digits(bytes, rug::integer::Order::MsfBe))
    }

    fn is_negative(&self) -> bool {
        self < &0
    }
}

#[cfg(feature = "openssl")]
#[cfg_attr(docsrs, doc(cfg(feature = "openssl")))]
impl Integer for openssl::bn::BigNum {
    fn to_be_bytes(&self) -> Vec<u8> {
        self.to_vec()
    }

    fn from_be_bytes(bytes: &[u8]) -> Option<Self> {
        Self::from_slice(bytes).ok()
    }

    fn is_negative(&self) -> bool {
        openssl::bn::BigNumRef::is_negative(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::string::ToString;

    #[test]
    fn natural_is_canonical_and_converts_to_backends() {
        let natural = Natural::from_be_bytes(&[0, 0, 1, 0, 1]);
        assert_eq!(natural.as_be_bytes(), &[1, 0, 1]);
        assert_eq!(natural.to_integer::<u64>(), Some(65_537));
        assert_eq!(natural.to_string(), "65537");
        assert_eq!(alloc::format!("{natural:?}"), "0x010001");
    }

    #[test]
    fn fixed_width_conversion_rejects_overflow() {
        assert_eq!(<u64 as Integer>::from_be_bytes(&[1; 9]), None);
        assert_eq!(<u64 as Integer>::from_be_bytes(&[1]), Some(1));
    }

    #[test]
    fn signed_backend_reports_negative_values() {
        let negative = BigInt::from(-1);
        assert!(Integer::is_negative(&negative));
        assert_eq!(Integer::to_be_bytes(&negative), alloc::vec![1]);
    }
}
