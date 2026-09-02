use alloc::vec::Vec;

use crate::{Error, Integer, Natural, Result};

/// A short Weierstrass curve `y² = x³ + ax + b (mod n)`.
#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Curve {
    /// Linear coefficient.
    pub a: Natural,
    /// Constant coefficient.
    pub b: Natural,
}

/// An affine elliptic-curve point.
#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Point {
    /// Affine x-coordinate.
    pub x: Natural,
    /// Affine y-coordinate.
    pub y: Natural,
}

/// One Atkin–Morain ECPP reduction from `n` to `q`.
#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EcppStep {
    /// Integer whose primality this step certifies, assuming `q` is prime.
    pub n: Natural,
    /// Curve over `Z/nZ`.
    pub curve: Curve,
    /// Point used by the elliptic Pocklington criterion.
    pub point: Point,
    /// Known curve order `m`.
    pub order: Natural,
    /// Prime divisor `q` of `m`, certified by the following proof node.
    pub q: Natural,
}

/// A node in a recursive primality certificate.
#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ProofNode {
    /// An ECPP reduction.
    EllipticCurve(EcppStep),
    /// A prime small enough for deterministic 64-bit Miller–Rabin.
    SmallPrime(u64),
}

/// A complete primality certificate ordered from the claimed prime to its base case.
#[derive(Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PrimalityProof {
    /// Certificate nodes. The last node must be [`ProofNode::SmallPrime`].
    pub nodes: Vec<ProofNode>,
}

impl PrimalityProof {
    /// Returns the integer certified by this proof.
    pub fn number(&self) -> Option<Natural> {
        match self.nodes.first()? {
            ProofNode::EllipticCurve(step) => Some(step.n.clone()),
            ProofNode::SmallPrime(prime) => Some(Natural::from_be_bytes(&prime.to_be_bytes())),
        }
    }

    /// Checks this proof without needing the original backend value.
    pub fn verify(&self) -> Result<()> {
        crate::prime::verify_proof(self)
    }

    /// Checks that this proof certifies `candidate`.
    pub fn verify_for<T: Integer>(&self, candidate: &T) -> Result<()> {
        if candidate.is_negative() {
            return Err(Error::InvalidInput("candidate must be non-negative"));
        }
        let expected = Natural::from_be_bytes(&candidate.to_be_bytes());
        if self.number().as_ref() != Some(&expected) {
            return Err(Error::InvalidProof("proof is for a different integer"));
        }
        self.verify()
    }
}

#[cfg(all(test, feature = "serde"))]
mod tests {
    use super::*;

    fn sample_proof() -> PrimalityProof {
        PrimalityProof {
            nodes: alloc::vec![ProofNode::SmallPrime(65_537)],
        }
    }

    #[test]
    fn postcard_round_trip() {
        let proof = sample_proof();
        let encoded = postcard::to_allocvec(&proof).unwrap();
        let decoded: PrimalityProof = postcard::from_bytes(&encoded).unwrap();
        assert_eq!(decoded, proof);
        decoded.verify().unwrap();
    }

    #[test]
    fn cbor_round_trip() {
        let proof = sample_proof();
        let encoded = serde_cbor_2::to_vec(&proof).unwrap();
        let decoded: PrimalityProof = serde_cbor_2::from_slice(&encoded).unwrap();
        assert_eq!(decoded, proof);
        decoded.verify().unwrap();
    }

    #[test]
    fn json_round_trip() {
        let proof = sample_proof();
        let encoded = serde_json::to_vec(&proof).unwrap();
        let decoded: PrimalityProof = serde_json::from_slice(&encoded).unwrap();
        assert_eq!(decoded, proof);
        decoded.verify().unwrap();
    }

    #[test]
    fn toml_round_trip() {
        let proof = sample_proof();
        let encoded = toml::to_string(&proof).unwrap();
        let decoded: PrimalityProof = toml::from_str(&encoded).unwrap();
        assert_eq!(decoded, proof);
        decoded.verify().unwrap();
    }

    #[test]
    fn yaml_round_trip() {
        let proof = sample_proof();
        let encoded = yaml_serde::to_string(&proof).unwrap();
        let decoded: PrimalityProof = yaml_serde::from_str(&encoded).unwrap();
        assert_eq!(decoded, proof);
        decoded.verify().unwrap();
    }
}
