#[derive(Clone, Copy, Debug)]
pub(crate) enum ClassPolynomial {
    // The embedded Hilbert class polynomials have degree at most two, so their
    // roots can be recovered directly instead of by general polynomial logic.
    Linear(i128),
    Quadratic { constant: i128, linear: i128 },
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Discriminant {
    pub(crate) value: i16,
    pub(crate) polynomial: ClassPolynomial,
}

/// Fundamental class-number-one and class-number-two discriminants excluding
/// the special `j = 0` and `j = 1728` cases.
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
