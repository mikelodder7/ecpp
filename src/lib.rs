#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]
#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

extern crate alloc;
#[cfg(feature = "std")]
extern crate std;

mod backend;
mod certificate;
mod cm;
mod error;
mod math;
pub mod prime;

pub use backend::{Integer, Natural};
pub use certificate::{Curve, EcppStep, Point, PrimalityProof, ProofNode};
pub use error::{Error, Result};
pub use prime::{ProvedPrime, ProverOptions};
