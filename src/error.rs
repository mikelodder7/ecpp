use crate::Natural;

/// Errors produced while constructing or checking a primality proof.
#[derive(Clone, Debug, Eq, PartialEq, thiserror::Error)]
pub enum Error {
    /// The input is negative, even, too small, or does not fit the requested backend.
    #[error("invalid input: {0}")]
    InvalidInput(&'static str),
    /// A compositeness witness was found.
    #[error("candidate is composite")]
    Composite,
    /// The configured CM/factor search did not find a proof step.
    #[error("ECPP search exhausted while proving {candidate}")]
    SearchExhausted {
        /// Candidate at which construction stopped.
        candidate: Natural,
    },
    /// The supplied certificate is malformed or fails a required identity.
    #[error("invalid primality proof: {0}")]
    InvalidProof(&'static str),
}

/// Result type used by this crate.
pub type Result<T> = core::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::string::ToString;

    #[test]
    fn errors_have_actionable_messages() {
        assert_eq!(Error::Composite.to_string(), "candidate is composite");
        assert_eq!(
            Error::InvalidInput("candidate must be non-negative").to_string(),
            "invalid input: candidate must be non-negative"
        );
        assert_eq!(
            Error::SearchExhausted {
                candidate: Natural::from_be_bytes(&[17]),
            }
            .to_string(),
            "ECPP search exhausted while proving 17"
        );
    }
}
