/// The largest class number in the embedded discriminant table.
pub(crate) const MAX_CLASS_NUMBER: usize = 8;

/// A signed integer stored as a big-endian magnitude, used for Hilbert class
/// polynomial coefficients too large for `i128`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SignedMagnitude {
    pub(crate) negative: bool,
    pub(crate) magnitude: &'static [u8],
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum ClassPolynomial {
    // Linear and quadratic Hilbert class polynomials are solved directly;
    // higher degrees carry their non-leading coefficients (constant term
    // first, monic leading coefficient implied) and are split with
    // randomized polynomial factoring.
    Linear(i128),
    Quadratic { constant: i128, linear: i128 },
    General(&'static [SignedMagnitude]),
}

impl ClassPolynomial {
    /// The class number of the discriminant, which equals the polynomial
    /// degree.
    pub(crate) const fn class_number(&self) -> u8 {
        match self {
            Self::Linear(_) => 1,
            Self::Quadratic { .. } => 2,
            Self::General(coefficients) => coefficients.len() as u8,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Discriminant {
    pub(crate) value: i16,
    pub(crate) polynomial: ClassPolynomial,
}

include!("cm_generated.rs");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_ordered_and_within_declared_class_numbers() {
        assert!(DISCRIMINANTS.len() > 300);
        let mut previous = 0u8;
        for discriminant in DISCRIMINANTS {
            let class_number = discriminant.polynomial.class_number();
            assert!((1..=8).contains(&class_number));
            assert!(class_number >= previous, "table ordered by class number");
            previous = class_number;
            assert!(discriminant.value < 0);
            if let ClassPolynomial::General(coefficients) = discriminant.polynomial {
                assert!(coefficients.len() >= 3);
                for coefficient in coefficients {
                    assert!(
                        !(coefficient.negative && coefficient.magnitude.is_empty()),
                        "negative zero coefficient"
                    );
                }
            }
        }
    }

    #[test]
    fn known_class_number_one_and_two_entries_are_present() {
        let linear = DISCRIMINANTS
            .iter()
            .find(|d| d.value == -163)
            .expect("D=-163 present");
        assert!(matches!(
            linear.polynomial,
            ClassPolynomial::Linear(-262_537_412_640_768_000)
        ));
        let quadratic = DISCRIMINANTS
            .iter()
            .find(|d| d.value == -427)
            .expect("D=-427 present");
        assert!(matches!(
            quadratic.polynomial,
            ClassPolynomial::Quadratic {
                constant: 155_041_756_222_618_916_546_936_832_000_000,
                linear: 15_611_455_512_523_783_919_812_608_000,
            }
        ));
    }
}
