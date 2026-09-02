#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Backend-neutral elliptic curve primality proving.
//!
//! `ecpp` constructs Atkin–Morain elliptic curve primality certificates and
//! verifies them independently of the integer type used by the caller. Unlike
//! a probable-prime test, a successful proof contains the data needed to check
//! the primality claim again without trusting the original proving run.
//!
//! # Proving and verifying
//!
//! ```
//! # #[cfg(all(feature = "num-bigint", feature = "getrandom"))]
//! # fn main() -> Result<(), ecpp::Error> {
//! use ecpp::prime;
//!
//! let candidate = 65_537u64;
//! let proof = prime::prove(&candidate)?;
//!
//! prime::verify(&candidate, &proof)?;
//! proof.verify()?;
//! assert_eq!(
//!     proof.number().and_then(|number| number.to_integer::<u64>()),
//!     Some(candidate),
//! );
//! # Ok(())
//! # }
//! # #[cfg(not(all(feature = "num-bigint", feature = "getrandom")))]
//! # fn main() {}
//! ```
//!
//! `prime::prove` accepts any type implementing `Integer`. The generic engine
//! uses an `ArithmeticBackend` selected independently from the input type.
//! Implementations are provided for `num-bigint`, fixed-width and heap-backed
//! `crypto-bigint`, `rug`, and OpenSSL. Other libraries can implement `Integer`
//! or use the canonical big-endian byte APIs.
//!
//! # Generating proved primes
//!
//! `prime::new` returns both a generated value and its `PrimalityProof`:
//!
//! ```
//! # #[cfg(all(feature = "num-bigint", feature = "getrandom"))]
//! # fn main() -> Result<(), ecpp::Error> {
//! use ecpp::prime;
//!
//! let generated = prime::new::<u128>(96)?;
//! prime::verify(&generated.prime, &generated.proof)?;
//! # Ok(())
//! # }
//! # #[cfg(not(all(feature = "num-bigint", feature = "getrandom")))]
//! # fn main() {}
//! ```
//!
//! Proving is a bounded search. `Error::SearchExhausted` means the configured
//! search bounds were reached; it is not evidence that the candidate is
//! composite. Use `prime::prove_with_options` to choose explicit bounds and
//! `prime::prove_with_rng` when the caller must supply the RNG.
//!
//! # Certificates and serialization
//!
//! Certificates contain backend-neutral `Natural` values. With the `serde`
//! feature, `PrimalityProof` can be serialized by any Serde data format.
//! Certificate integers are lowercase hexadecimal strings in human-readable
//! formats and byte strings in binary formats. A deserialized proof must still
//! be checked with `PrimalityProof::verify` or `PrimalityProof::verify_for`.
//!
//! Heap-backed proving is `no_std` with the `num-bigint` or `rug` engine. The
//! `fixed` module instead uses `crypto_bigint::Uint` throughout and stores
//! proofs in a caller-selected number of stack slots, so it works with
//! `--no-default-features --features crypto-bigint` and no allocator. With
//! `alloc`, `fixed::PrimalityProof::to_alloc` converts that stack certificate
//! into the same neutral format verified by every heap engine.

#[cfg(feature = "alloc")]
extern crate alloc;
#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "alloc")]
pub mod arithmetic;
#[cfg(feature = "alloc")]
mod backend;
#[cfg(feature = "alloc")]
mod certificate;
#[cfg(any(feature = "alloc", feature = "crypto-bigint"))]
mod cm;
#[cfg(feature = "alloc")]
mod engine;
#[cfg(feature = "alloc")]
mod error;
#[cfg(feature = "crypto-bigint")]
pub mod fixed;
#[cfg(feature = "alloc")]
pub mod prime;
mod proved_prime;

#[cfg(feature = "alloc")]
pub use backend::{Integer, Natural};
#[cfg(feature = "alloc")]
pub use certificate::{Curve, EcppStep, Point, PrimalityProof, ProofNode};
#[cfg(feature = "alloc")]
pub use error::{Error, Result};
#[cfg(feature = "crypto-bigint")]
pub use fixed::PrimalityProof as FixedPrimalityProof;
#[cfg(feature = "alloc")]
pub use prime::ProverOptions;
pub use proved_prime::ProvedPrime;
