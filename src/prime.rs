//! Prime generation, proving, and certificate verification.

use alloc::vec;

use rand_core::CryptoRng;

use crate::arithmetic::ArithmeticBackend;
use crate::{Error, Integer, PrimalityProof, ProvedPrime, Result};

pub use crate::engine::ProverOptions;

/// Generates a proved prime with a selected arithmetic backend and RNG.
pub fn from_rng_with_backend<B: ArithmeticBackend, T: Integer, R: CryptoRng + ?Sized>(
    bit_length: usize,
    rng: &mut R,
) -> Result<ProvedPrime<T, PrimalityProof>> {
    if bit_length < 2 {
        return Err(Error::InvalidInput("bit length must be at least two"));
    }
    let byte_length = bit_length.div_ceil(8);
    let excess = byte_length * 8 - bit_length;
    loop {
        let mut bytes = vec![0u8; byte_length];
        rng.fill_bytes(&mut bytes);
        bytes[0] &= u8::MAX >> excess;
        bytes[0] |= 1u8 << (7 - excess);
        bytes[byte_length - 1] |= 1;
        let proof = match crate::engine::prove::<B, R>(&bytes, rng, ProverOptions::default()) {
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

/// Generates a proved prime using the default enabled backend and `rng`.
#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
pub fn from_rng<T: Integer, R: CryptoRng + ?Sized>(
    bit_length: usize,
    rng: &mut R,
) -> Result<ProvedPrime<T, PrimalityProof>> {
    #[cfg(feature = "num-bigint")]
    return from_rng_with_backend::<crate::arithmetic::NumBigint, T, R>(bit_length, rng);
    #[cfg(all(not(feature = "num-bigint"), feature = "crypto-bigint"))]
    return from_rng_with_backend::<crate::arithmetic::CryptoBigint, T, R>(bit_length, rng);
    #[cfg(all(
        not(feature = "num-bigint"),
        not(feature = "crypto-bigint"),
        feature = "rug"
    ))]
    return from_rng_with_backend::<crate::arithmetic::Rug, T, R>(bit_length, rng);
    #[cfg(all(
        not(feature = "num-bigint"),
        not(feature = "crypto-bigint"),
        not(feature = "rug"),
        feature = "openssl"
    ))]
    from_rng_with_backend::<crate::arithmetic::OpenSsl, T, R>(bit_length, rng)
}

/// Generates a proved prime using the default backend and operating-system RNG.
#[cfg(all(
    feature = "getrandom",
    any(
        feature = "num-bigint",
        feature = "crypto-bigint",
        feature = "rug",
        feature = "openssl"
    )
))]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn new<T: Integer>(bit_length: usize) -> Result<ProvedPrime<T, PrimalityProof>> {
    from_rng(bit_length, &mut rand_core::UnwrapErr(getrandom::SysRng))
}

/// Constructs an ECPP proof with a selected backend, RNG, and search bounds.
pub fn prove_with_backend<B: ArithmeticBackend, T: Integer, R: CryptoRng + ?Sized>(
    candidate: &T,
    rng: &mut R,
    options: ProverOptions,
) -> Result<PrimalityProof> {
    if candidate.is_negative() {
        return Err(Error::InvalidInput("candidate must be non-negative"));
    }
    crate::engine::prove::<B, R>(&candidate.to_be_bytes(), rng, options)
}

/// Constructs an ECPP proof with explicit search bounds and the default backend.
#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
pub fn prove_with_options<T: Integer, R: CryptoRng + ?Sized>(
    candidate: &T,
    rng: &mut R,
    options: ProverOptions,
) -> Result<PrimalityProof> {
    #[cfg(feature = "num-bigint")]
    return prove_with_backend::<crate::arithmetic::NumBigint, T, R>(candidate, rng, options);
    #[cfg(all(not(feature = "num-bigint"), feature = "crypto-bigint"))]
    return prove_with_backend::<crate::arithmetic::CryptoBigint, T, R>(candidate, rng, options);
    #[cfg(all(
        not(feature = "num-bigint"),
        not(feature = "crypto-bigint"),
        feature = "rug"
    ))]
    return prove_with_backend::<crate::arithmetic::Rug, T, R>(candidate, rng, options);
    #[cfg(all(
        not(feature = "num-bigint"),
        not(feature = "crypto-bigint"),
        not(feature = "rug"),
        feature = "openssl"
    ))]
    prove_with_backend::<crate::arithmetic::OpenSsl, T, R>(candidate, rng, options)
}

/// Constructs an ECPP proof using `rng` and the default backend.
#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
pub fn prove_with_rng<T: Integer, R: CryptoRng + ?Sized>(
    candidate: &T,
    rng: &mut R,
) -> Result<PrimalityProof> {
    prove_with_options(candidate, rng, ProverOptions::default())
}

/// Constructs an ECPP proof with the default backend and operating-system RNG.
#[cfg(all(
    feature = "getrandom",
    any(
        feature = "num-bigint",
        feature = "crypto-bigint",
        feature = "rug",
        feature = "openssl"
    )
))]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn prove<T: Integer>(candidate: &T) -> Result<PrimalityProof> {
    prove_with_rng(candidate, &mut rand_core::UnwrapErr(getrandom::SysRng))
}

/// Constructs a proof from an unsigned magnitude with a selected backend.
pub fn prove_be_bytes_with_backend<B: ArithmeticBackend, R: CryptoRng + ?Sized>(
    candidate: &[u8],
    rng: &mut R,
    options: ProverOptions,
) -> Result<PrimalityProof> {
    crate::engine::prove::<B, R>(candidate, rng, options)
}

/// Constructs a proof from an unsigned magnitude using the default backend.
#[cfg(all(
    feature = "getrandom",
    any(
        feature = "num-bigint",
        feature = "crypto-bigint",
        feature = "rug",
        feature = "openssl"
    )
))]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn prove_be_bytes(candidate: &[u8]) -> Result<PrimalityProof> {
    prove(&crate::Natural::from_be_bytes(candidate))
}

/// Verifies a proof against `candidate` with a selected backend.
pub fn verify_with_backend<B: ArithmeticBackend, T: Integer>(
    candidate: &T,
    proof: &PrimalityProof,
) -> Result<()> {
    proof.verify_for_with::<B, T>(candidate)
}

/// Verifies a proof against `candidate` with the default backend.
#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
pub fn verify<T: Integer>(candidate: &T, proof: &PrimalityProof) -> Result<()> {
    proof.verify_for(candidate)
}

/// Verifies a proof against an unsigned magnitude with a selected backend.
pub fn verify_be_bytes_with_backend<B: ArithmeticBackend>(
    candidate: &[u8],
    proof: &PrimalityProof,
) -> Result<()> {
    let expected = crate::Natural::from_be_bytes(candidate);
    if proof.number().as_ref() != Some(&expected) {
        return Err(Error::InvalidProof("proof is for a different integer"));
    }
    proof.verify_with::<B>()
}

/// Verifies a proof against an unsigned magnitude with the default backend.
#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
pub fn verify_be_bytes(candidate: &[u8], proof: &PrimalityProof) -> Result<()> {
    verify(&crate::Natural::from_be_bytes(candidate), proof)
}

/// Proves primality and returns `false` when construction fails.
#[cfg(all(
    feature = "getrandom",
    any(
        feature = "num-bigint",
        feature = "crypto-bigint",
        feature = "rug",
        feature = "openssl"
    )
))]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn check<T: Integer>(candidate: &T) -> bool {
    prove(candidate).is_ok()
}

/// RNG-explicit counterpart to `check`.
#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
pub fn check_with<T: Integer, R: CryptoRng + ?Sized>(candidate: &T, rng: &mut R) -> bool {
    prove_with_rng(candidate, rng).is_ok()
}

/// Alias for `check`; a successful ECPP proof is stronger than screening.
#[cfg(all(
    feature = "getrandom",
    any(
        feature = "num-bigint",
        feature = "crypto-bigint",
        feature = "rug",
        feature = "openssl"
    )
))]
#[cfg_attr(docsrs, doc(cfg(feature = "getrandom")))]
pub fn strong_check<T: Integer>(candidate: &T) -> bool {
    check(candidate)
}

/// RNG-explicit counterpart to `strong_check`.
#[cfg(any(
    feature = "num-bigint",
    feature = "crypto-bigint",
    feature = "rug",
    feature = "openssl"
))]
pub fn strong_check_with<T: Integer, R: CryptoRng + ?Sized>(candidate: &T, rng: &mut R) -> bool {
    check_with(candidate, rng)
}

#[cfg(all(
    test,
    any(
        feature = "num-bigint",
        feature = "crypto-bigint",
        feature = "rug",
        feature = "openssl"
    )
))]
mod tests {
    use super::*;
    use rand::{SeedableRng, rngs::StdRng};

    #[test]
    #[cfg(feature = "num-bigint")]
    fn proves_and_verifies_a_multistep_prime() {
        let prime =
            num_bigint::BigUint::parse_bytes(b"340282366920938463463374607431768211297", 10)
                .unwrap();
        let mut rng = StdRng::seed_from_u64(7);
        let proof = prove_with_rng(&prime, &mut rng).unwrap();
        assert!(proof.nodes.len() >= 2);
        proof.verify_for(&prime).unwrap();

        #[cfg(feature = "rug")]
        proof.verify_with::<crate::arithmetic::Rug>().unwrap();
        #[cfg(feature = "openssl")]
        proof.verify_with::<crate::arithmetic::OpenSsl>().unwrap();

        let mut changed = proof;
        if let crate::ProofNode::EllipticCurve(step) = &mut changed.nodes[0] {
            step.point.x = crate::Natural::from_be_bytes(&[1]);
        }
        assert!(changed.verify().is_err());
    }

    #[test]
    #[ignore = "expensive 256-bit ECPP regression"]
    #[cfg(feature = "num-bigint")]
    fn proves_secp256k1_field_modulus() {
        let prime = num_bigint::BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap();
        let mut rng = StdRng::seed_from_u64(29);
        let proof = prove_with_rng(&prime, &mut rng).unwrap();
        proof.verify_for(&prime).unwrap();
    }

    #[test]
    fn rejects_a_carmichael_number() {
        let mut rng = StdRng::seed_from_u64(11);
        assert_eq!(prove_with_rng(&561u64, &mut rng), Err(Error::Composite));
    }

    #[test]
    fn generates_a_proved_prime() {
        let mut rng = StdRng::seed_from_u64(19);
        let generated: ProvedPrime<u64, PrimalityProof> = from_rng(16, &mut rng).unwrap();
        generated.proof.verify_for(&generated.prime).unwrap();
        assert!(check_with(&65_537u64, &mut rng));
        assert!(strong_check_with(&65_537u64, &mut rng));
        assert!(!check_with(&561u64, &mut rng));
    }

    #[test]
    #[cfg(feature = "getrandom")]
    fn operating_system_rng_apis() {
        let generated: ProvedPrime<u64, PrimalityProof> = new(16).unwrap();
        generated.proof.verify_for(&generated.prime).unwrap();
        let proof = prove(&65_537u64).unwrap();
        verify(&65_537u64, &proof).unwrap();
        assert!(check(&65_537u64));
        assert!(strong_check(&65_537u64));
    }

    #[test]
    #[cfg(feature = "num-bigint")]
    fn explicit_backend_and_byte_apis() {
        let candidate = 65_537u64;
        let bytes = candidate.to_be_bytes();
        let mut rng = StdRng::seed_from_u64(23);
        let proof = prove_be_bytes_with_backend::<crate::arithmetic::NumBigint, _>(
            &bytes,
            &mut rng,
            ProverOptions::default(),
        )
        .unwrap();
        verify_with_backend::<crate::arithmetic::NumBigint, _>(&candidate, &proof).unwrap();
        verify_be_bytes_with_backend::<crate::arithmetic::NumBigint>(&bytes, &proof).unwrap();

        let generated: ProvedPrime<u64, PrimalityProof> =
            from_rng_with_backend::<crate::arithmetic::NumBigint, _, _>(16, &mut rng).unwrap();
        generated
            .proof
            .verify_for_with::<crate::arithmetic::NumBigint, _>(&generated.prime)
            .unwrap();
    }

    #[test]
    #[cfg(feature = "crypto-bigint")]
    fn crypto_bigint_backend_proves_without_num_bigint() {
        let value = 340_282_366_920_938_463_463_374_607_431_768_211_297u128;
        let bytes = value.to_be_bytes();
        let candidate = crypto_bigint::U256::from_u128(value);
        let boxed_candidate = crypto_bigint::BoxedUint::from_be_slice_vartime(&bytes);
        let mut rng = StdRng::seed_from_u64(7);
        let proof = prove_with_backend::<crate::arithmetic::CryptoBigint, _, _>(
            &candidate,
            &mut rng,
            ProverOptions::default(),
        )
        .unwrap();
        proof
            .verify_for_with::<crate::arithmetic::CryptoBigint, _>(&candidate)
            .unwrap();
        proof
            .verify_for_with::<crate::arithmetic::CryptoBigint, _>(&boxed_candidate)
            .unwrap();

        let uint_proof = prove_with_backend::<
            crate::arithmetic::CryptoUint<{ crypto_bigint::U512::LIMBS }>,
            _,
            _,
        >(&candidate, &mut rng, ProverOptions::default())
        .unwrap();
        uint_proof
            .verify_for_with::<crate::arithmetic::CryptoUint<{ crypto_bigint::U512::LIMBS }>, _>(
                &candidate,
            )
            .unwrap();
    }

    #[test]
    #[cfg(feature = "rug")]
    fn rug_backend_proves_without_num_bigint() {
        let candidate = rug::Integer::from(340_282_366_920_938_463_463_374_607_431_768_211_297u128);
        let mut rng = StdRng::seed_from_u64(7);
        let proof = prove_with_backend::<crate::arithmetic::Rug, _, _>(
            &candidate,
            &mut rng,
            ProverOptions::default(),
        )
        .unwrap();
        proof.verify_for(&candidate).unwrap();
    }

    #[test]
    #[cfg(feature = "openssl")]
    fn openssl_backend_proves_without_num_bigint() {
        let bytes = 340_282_366_920_938_463_463_374_607_431_768_211_297u128.to_be_bytes();
        let candidate = openssl::bn::BigNum::from_slice(&bytes).unwrap();
        let mut rng = StdRng::seed_from_u64(7);
        let proof = prove_with_backend::<crate::arithmetic::OpenSsl, _, _>(
            &candidate,
            &mut rng,
            ProverOptions::default(),
        )
        .unwrap();
        proof.verify_for(&candidate).unwrap();
    }
}
