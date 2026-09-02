/// A generated prime together with the proof that certifies it.
///
/// `T` is the caller's integer backend and `P` is the independently chosen
/// certificate representation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProvedPrime<T, P> {
    /// Generated prime value.
    pub prime: T,
    /// Certificate proving that `prime` is prime.
    pub proof: P,
}

impl<T, P> ProvedPrime<T, P> {
    /// Splits the generated value into its prime and proof.
    pub fn into_parts(self) -> (T, P) {
        (self.prime, self.proof)
    }
}
