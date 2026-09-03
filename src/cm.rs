#[derive(Clone, Copy, Debug)]
pub(crate) enum ClassPolynomial {
    // The embedded Hilbert class polynomials have degree at most three.
    // Linear and quadratic roots are recovered directly; cubic roots use
    // randomized polynomial splitting.
    Linear(i128),
    Quadratic {
        constant: i128,
        linear: i128,
    },
    Cubic {
        constant: i128,
        linear: i128,
        quadratic: i128,
    },
}

impl ClassPolynomial {
    /// The class number of the discriminant, which equals the polynomial
    /// degree.
    pub(crate) const fn class_number(&self) -> u8 {
        match self {
            Self::Linear(_) => 1,
            Self::Quadratic { .. } => 2,
            Self::Cubic { .. } => 3,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Discriminant {
    pub(crate) value: i16,
    pub(crate) polynomial: ClassPolynomial,
}

/// Fundamental discriminants of class number one, two, and three excluding
/// the special `j = 0` and `j = 1728` cases.
///
/// The class-number-one and class-number-two lists are complete. The
/// class-number-three list contains the twelve of sixteen fundamental
/// discriminants whose Hilbert class polynomial coefficients fit `i128`;
/// `-499`, `-643`, `-883`, and `-907` are omitted because theirs do not.
pub(crate) const DISCRIMINANTS: [Discriminant; 37] = [
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
    Discriminant {
        value: -23,
        polynomial: ClassPolynomial::Cubic {
            constant: 12_771_880_859_375,
            linear: -5_151_296_875,
            quadratic: 3_491_750,
        },
    },
    Discriminant {
        value: -31,
        polynomial: ClassPolynomial::Cubic {
            constant: 1_566_028_350_940_383,
            linear: -58_682_638_134,
            quadratic: 39_491_307,
        },
    },
    Discriminant {
        value: -59,
        polynomial: ClassPolynomial::Cubic {
            constant: 374_643_194_001_883_136,
            linear: -140_811_576_541_184,
            quadratic: 30_197_678_080,
        },
    },
    Discriminant {
        value: -83,
        polynomial: ClassPolynomial::Cubic {
            constant: 549_755_813_888_000_000_000,
            linear: -41_490_055_168_000_000,
            quadratic: 2_691_907_584_000,
        },
    },
    Discriminant {
        value: -107,
        polynomial: ClassPolynomial::Cubic {
            constant: 337_618_789_203_968_000_000_000,
            linear: -6_764_523_159_552_000_000,
            quadratic: 129_783_279_616_000,
        },
    },
    Discriminant {
        value: -139,
        polynomial: ClassPolynomial::Cubic {
            constant: 67_408_489_017_571_610_198_016,
            linear: -53_041_786_755_137_667_072,
            quadratic: 12_183_160_834_031_616,
        },
    },
    Discriminant {
        value: -211,
        polynomial: ClassPolynomial::Cubic {
            constant: 5_310_823_021_408_898_698_117_644_288,
            linear: 277_390_576_406_111_100_862_464,
            quadratic: 65_873_587_288_630_099_968,
        },
    },
    Discriminant {
        value: -283,
        polynomial: ClassPolynomial::Cubic {
            constant: 201_371_843_156_955_365_376_000_000_000,
            linear: 90_839_236_535_446_929_408_000_000,
            quadratic: 89_611_323_386_832_801_792_000,
        },
    },
    Discriminant {
        value: -307,
        polynomial: ClassPolynomial::Cubic {
            constant: 8_987_619_631_060_626_702_336_000_000_000,
            linear: -5_083_646_425_734_146_162_688_000_000,
            quadratic: 805_016_812_009_981_390_848_000,
        },
    },
    Discriminant {
        value: -331,
        polynomial: ClassPolynomial::Cubic {
            constant: 56_176_242_840_389_398_230_218_488_594_563_072,
            linear: 368_729_929_041_040_103_875_232_661_504,
            quadratic: 6_647_404_730_173_793_386_463_232,
        },
    },
    Discriminant {
        value: -379,
        polynomial: ClassPolynomial::Cubic {
            constant: 15_443_600_047_689_011_948_024_601_807_415_148_544,
            linear: -121_567_791_009_880_876_719_538_528_321_536,
            quadratic: 364_395_404_104_624_239_018_246_144,
        },
    },
    Discriminant {
        value: -547,
        polynomial: ClassPolynomial::Cubic {
            constant: 83_303_937_570_678_403_968_635_240_448_000_000_000,
            linear: -139_712_328_431_787_827_943_469_744_128_000_000,
            quadratic: 81_297_395_539_631_654_721_637_478_400_000,
        },
    },
];
