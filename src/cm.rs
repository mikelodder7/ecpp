use alloc::vec::Vec;

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::math::{integer_sqrt, is_square, mod_signed, mod_sub, modular_sqrt};

#[derive(Clone, Copy, Debug)]
pub(crate) enum ClassPolynomial {
    Linear(i128),
    Quadratic { constant: i128, linear: i128 },
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Discriminant {
    pub(crate) value: i16,
    pub(crate) polynomial: ClassPolynomial,
}

/// Fundamental class-number-one and class-number-two discriminants excluding
/// the special j=0 and j=1728 cases.
pub(crate) const DISCRIMINANTS: [Discriminant; 25] = [
    Discriminant {
        value: -7,
        polynomial: ClassPolynomial::Linear(-3_375),
    },
    Discriminant {
        value: -8,
        polynomial: ClassPolynomial::Linear(8_000),
    },
    Discriminant {
        value: -11,
        polynomial: ClassPolynomial::Linear(-32_768),
    },
    Discriminant {
        value: -19,
        polynomial: ClassPolynomial::Linear(-884_736),
    },
    Discriminant {
        value: -43,
        polynomial: ClassPolynomial::Linear(-884_736_000),
    },
    Discriminant {
        value: -67,
        polynomial: ClassPolynomial::Linear(-147_197_952_000),
    },
    Discriminant {
        value: -163,
        polynomial: ClassPolynomial::Linear(-262_537_412_640_768_000),
    },
    Discriminant {
        value: -15,
        polynomial: ClassPolynomial::Quadratic {
            constant: -121_287_375,
            linear: 191_025,
        },
    },
    Discriminant {
        value: -20,
        polynomial: ClassPolynomial::Quadratic {
            constant: -681_472_000,
            linear: -1_264_000,
        },
    },
    Discriminant {
        value: -24,
        polynomial: ClassPolynomial::Quadratic {
            constant: 14_670_139_392,
            linear: -4_834_944,
        },
    },
    Discriminant {
        value: -35,
        polynomial: ClassPolynomial::Quadratic {
            constant: -134_217_728_000,
            linear: 117_964_800,
        },
    },
    Discriminant {
        value: -40,
        polynomial: ClassPolynomial::Quadratic {
            constant: 9_103_145_472_000,
            linear: -425_692_800,
        },
    },
    Discriminant {
        value: -51,
        polynomial: ClassPolynomial::Quadratic {
            constant: 6_262_062_317_568,
            linear: 5_541_101_568,
        },
    },
    Discriminant {
        value: -52,
        polynomial: ClassPolynomial::Quadratic {
            constant: -567_663_552_000_000,
            linear: -6_896_880_000,
        },
    },
    Discriminant {
        value: -88,
        polynomial: ClassPolynomial::Quadratic {
            constant: 15_798_135_578_688_000_000,
            linear: -6_294_842_640_000,
        },
    },
    Discriminant {
        value: -91,
        polynomial: ClassPolynomial::Quadratic {
            constant: -3_845_689_020_776_448,
            linear: 10_359_073_013_760,
        },
    },
    Discriminant {
        value: -115,
        polynomial: ClassPolynomial::Quadratic {
            constant: 130_231_327_260_672_000,
            linear: 427_864_611_225_600,
        },
    },
    Discriminant {
        value: -123,
        polynomial: ClassPolynomial::Quadratic {
            constant: 148_809_594_175_488_000_000,
            linear: 1_354_146_840_576_000,
        },
    },
    Discriminant {
        value: -148,
        polynomial: ClassPolynomial::Quadratic {
            constant: -7_898_242_515_936_467_904_000_000,
            linear: -39_660_183_801_072_000,
        },
    },
    Discriminant {
        value: -187,
        polynomial: ClassPolynomial::Quadratic {
            constant: -3_845_689_020_776_448_000_000,
            linear: 4_545_336_381_788_160_000,
        },
    },
    Discriminant {
        value: -232,
        polynomial: ClassPolynomial::Quadratic {
            constant: 14_871_070_713_157_137_145_512_000_000_000,
            linear: -604_729_957_849_891_344_000,
        },
    },
    Discriminant {
        value: -235,
        polynomial: ClassPolynomial::Quadratic {
            constant: 11_946_621_170_462_723_407_872_000,
            linear: 823_177_419_449_425_920_000,
        },
    },
    Discriminant {
        value: -267,
        polynomial: ClassPolynomial::Quadratic {
            constant: 531_429_662_672_621_376_897_024_000_000,
            linear: 19_683_091_854_079_488_000_000,
        },
    },
    Discriminant {
        value: -403,
        polynomial: ClassPolynomial::Quadratic {
            constant: -108_844_203_402_491_055_833_088_000_000,
            linear: 2_452_811_389_229_331_391_979_520_000,
        },
    },
    Discriminant {
        value: -427,
        polynomial: ClassPolynomial::Quadratic {
            constant: 155_041_756_222_618_916_546_936_832_000_000,
            linear: 15_611_455_512_523_783_919_812_608_000,
        },
    },
];

pub(crate) fn j_invariants(
    candidate: &BigUint,
    polynomial: ClassPolynomial,
) -> Option<Vec<BigUint>> {
    match polynomial {
        ClassPolynomial::Linear(root) => Some(alloc::vec![mod_signed(root, candidate)]),
        ClassPolynomial::Quadratic { constant, linear } => {
            let linear_mod = mod_signed(linear, candidate);
            let constant_mod = mod_signed(constant, candidate);
            let discriminant = mod_sub(
                &(&linear_mod * &linear_mod % candidate),
                &(&constant_mod * 4u8 % candidate),
                candidate,
            );
            let square_root = modular_sqrt(&discriminant, candidate)?;
            let inverse_two = (candidate + 1u8) >> 1usize;
            let minus_linear = if linear_mod.is_zero() {
                BigUint::zero()
            } else {
                candidate - linear_mod
            };
            let first =
                (mod_sub(&minus_linear, &square_root, candidate) * &inverse_two) % candidate;
            let second = ((minus_linear + square_root) * inverse_two) % candidate;
            Some(alloc::vec![first, second])
        }
    }
}

/// Cohen's modified Cornacchia algorithm for `4n = u² + |D|v²`.
pub(crate) fn cornacchia(candidate: &BigUint, discriminant: i16) -> Option<(BigUint, BigUint)> {
    let absolute = discriminant.unsigned_abs();
    let residue = mod_signed(discriminant as i128, candidate);
    let mut root = modular_sqrt(&residue, candidate)?;
    let expected_odd = absolute & 1 == 1;
    let root_odd = (&root & BigUint::one()) == BigUint::one();
    if root_odd != expected_odd {
        root = candidate - root;
    }

    let mut previous = candidate << 1usize;
    let mut current = root;
    let limit = integer_sqrt(&(candidate << 2usize));
    while current > limit {
        let remainder = &previous % &current;
        previous = current;
        current = remainder;
        if current.is_zero() {
            return None;
        }
    }

    let four_n = candidate << 2usize;
    let square = &current * &current;
    if square > four_n {
        return None;
    }
    let remainder = four_n - square;
    if &remainder % absolute != BigUint::zero() {
        return None;
    }
    let v_squared = remainder / absolute;
    let v = is_square(&v_squared)?;
    Some((current, v))
}
