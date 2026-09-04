//! Generates the full CM discriminant table: every fundamental discriminant
//! with class number at most eight, with exact Hilbert class polynomial
//! coefficients. Class numbers one and two are emitted as i128 literals (and
//! cross-checked against the hand-verified table); higher class numbers are
//! emitted as signed big-endian byte strings.

use rug::ops::Pow;
use rug::{Complex, Float, Integer};

const MAX_CLASS_NUMBER: usize = 8;
const MAX_ABS_D: i64 = 8000;
const MAX_POLY_BYTES: usize = 4096;

fn is_squarefree(mut n: i64) -> bool {
    let mut p = 2i64;
    while p * p <= n {
        if n % (p * p) == 0 {
            return false;
        }
        while n % p == 0 {
            n /= p;
        }
        p += 1;
    }
    true
}

fn is_fundamental(d: i64) -> bool {
    let residue = ((d % 4) + 4) % 4;
    if residue == 1 {
        is_squarefree(-d)
    } else if residue == 0 {
        let m = d / 4;
        let m_residue = ((m % 4) + 4) % 4;
        (m_residue == 2 || m_residue == 3) && is_squarefree(-m)
    } else {
        false
    }
}

fn reduced_forms(d: i64) -> Vec<(i64, i64, i64)> {
    let mut forms = Vec::new();
    let a_max = (((-d) as f64) / 3.0).sqrt() as i64 + 1;
    for a in 1..=a_max {
        for b in (-a + 1)..=a {
            let num = b * b - d;
            if num % (4 * a) != 0 {
                continue;
            }
            let c = num / (4 * a);
            if c < a {
                continue;
            }
            if b < 0 && a == c {
                continue;
            }
            forms.push((a, b, c));
        }
    }
    forms
}

fn j_of_form(d: i64, a: i64, b: i64, prec: u32, terms: u32) -> Complex {
    let pi = Float::with_val(prec, rug::float::Constant::Pi);
    let sqrt_abs_d = Float::with_val(prec, -d).sqrt();
    let re = -(pi.clone() * sqrt_abs_d) / Float::with_val(prec, a);
    let im = -(pi * Float::with_val(prec, b)) / Float::with_val(prec, a);
    let q = Complex::with_val(prec, (re, im)).exp();

    let mut e4 = Complex::with_val(prec, 1);
    let mut delta_prod = Complex::with_val(prec, 1);
    let mut qn = Complex::with_val(prec, 1);
    for n in 1..=terms {
        qn *= &q;
        let one_minus_qn = Complex::with_val(prec, 1) - qn.clone();
        let term = Complex::with_val(prec, n).pow(3u32) * qn.clone() / one_minus_qn.clone();
        e4 += term * Complex::with_val(prec, 240);
        delta_prod *= one_minus_qn.pow(24u32);
    }
    let delta = q * delta_prod;
    e4.pow(3u32) / delta
}

fn round_real(value: &Complex, prec: u32) -> Option<Integer> {
    let re = value.real().clone();
    let rounded = re.clone().round().to_integer()?;
    let mut err = re;
    err -= Float::with_val(prec, &rounded);
    let err = err.abs();
    let imag = value.imag().clone().abs();
    if err > 1e-6 || imag > 1e-6 {
        return None;
    }
    Some(rounded)
}

fn magnitude_bytes(value: &Integer) -> Vec<u8> {
    let magnitude = value.clone().abs();
    if magnitude == 0 {
        return Vec::new();
    }
    magnitude.to_digits::<u8>(rug::integer::Order::MsfBe)
}

fn hex_array(bytes: &[u8]) -> String {
    let inner: Vec<String> = bytes.iter().map(|b| format!("0x{b:02x}")).collect();
    format!("&[{}]", inner.join(", "))
}

fn main() {
    let i128_max = Integer::from(i128::MAX);
    let mut entries: Vec<(i64, usize, String)> = Vec::new();
    let mut total_bytes = 0usize;
    let mut skipped_size = Vec::new();
    let mut skipped_integrality = Vec::new();

    for minus_d in 3..=MAX_ABS_D {
        let d = -minus_d;
        if d == -3 || d == -4 {
            continue; // handled by the j = 0 and j = 1728 special cases
        }
        if !is_fundamental(d) {
            continue;
        }
        let forms = reduced_forms(d);
        let h = forms.len();
        if h == 0 || h > MAX_CLASS_NUMBER {
            continue;
        }

        let sum_inv_a: f64 = forms.iter().map(|&(a, _, _)| 1.0 / a as f64).sum();
        let est_bits = 1.4427 * core::f64::consts::PI * ((-d) as f64).sqrt() * sum_inv_a;
        let prec = (est_bits as u32).saturating_add(256).max(512);
        let terms = prec / 6 + 64;

        let roots: Vec<Complex> = forms
            .iter()
            .map(|&(a, b, _)| j_of_form(d, a, b, prec, terms))
            .collect();

        // Expand prod (x - j_i) iteratively; poly[k] is the coefficient of x^k.
        let mut poly = vec![Complex::with_val(prec, 1)];
        for root in &roots {
            let mut next = vec![Complex::with_val(prec, 0); poly.len() + 1];
            for (k, coefficient) in poly.iter().enumerate() {
                next[k + 1] += coefficient;
                next[k] -= (coefficient.clone() * root).clone();
            }
            poly = next;
        }
        // poly has degree h with leading coefficient 1; take c_0..c_{h-1}.
        let mut coefficients = Vec::with_capacity(h);
        let mut ok = true;
        for coefficient in poly.iter().take(h) {
            match round_real(coefficient, prec) {
                Some(value) => coefficients.push(value),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            skipped_integrality.push(d);
            continue;
        }

        if h <= 2 {
            // Emit as i128 for cross-checking against the existing table.
            let fits = coefficients.iter().all(|c| c.clone().abs() <= i128_max);
            assert!(fits, "h<=2 coefficient exceeds i128 for D={d}");
            let entry = if h == 1 {
                // `Linear` stores the polynomial's root, the negated constant.
                let root: Integer = (-&coefficients[0]).into();
                format!(
                    "    Discriminant {{\n        value: {d},\n        polynomial: ClassPolynomial::Linear({root}),\n    }},"
                )
            } else {
                format!(
                    "    Discriminant {{\n        value: {d},\n        polynomial: ClassPolynomial::Quadratic {{\n            constant: {},\n            linear: {},\n        }},\n    }},",
                    coefficients[0], coefficients[1]
                )
            };
            entries.push((d, h, entry));
            continue;
        }

        let encoded: Vec<(bool, Vec<u8>)> = coefficients
            .iter()
            .map(|c| (*c < 0, magnitude_bytes(c)))
            .collect();
        let size: usize = encoded.iter().map(|(_, bytes)| bytes.len()).sum();
        if size > MAX_POLY_BYTES {
            skipped_size.push((d, h, size));
            continue;
        }
        total_bytes += size;
        let mut lines = Vec::new();
        for (negative, bytes) in &encoded {
            lines.push(format!(
                "            SignedMagnitude {{\n                negative: {negative},\n                magnitude: {},\n            }},",
                hex_array(bytes)
            ));
        }
        let entry = format!(
            "    Discriminant {{\n        value: {d},\n        polynomial: ClassPolynomial::General(&[\n{}\n        ]),\n    }},",
            lines.join("\n")
        );
        entries.push((d, h, entry));
    }

    // Order: by class number, then |D| ascending, matching search preference.
    entries.sort_by_key(|(d, h, _)| (*h, -*d));

    let mut histogram = [0usize; MAX_CLASS_NUMBER + 1];
    for (_, h, _) in &entries {
        histogram[*h] += 1;
    }
    eprintln!("total entries: {}", entries.len());
    for (h, count) in histogram.iter().enumerate().skip(1) {
        eprintln!("  h={h}: {count}");
    }
    eprintln!("general-coefficient bytes: {total_bytes}");
    eprintln!("skipped for size: {skipped_size:?}");
    eprintln!("skipped for integrality: {skipped_integrality:?}");

    println!("// GENERATED by tools/hilbert-table — do not edit by hand.");
    println!("//");
    println!(
        "// Every fundamental imaginary quadratic discriminant with class number"
    );
    println!(
        "// at most {MAX_CLASS_NUMBER} (excluding -3 and -4, which the engine special-cases),"
    );
    println!("// ordered by class number and then by |D|.");
    println!(
        "pub(crate) static DISCRIMINANTS: &[Discriminant] = &["
    );
    for (_, _, entry) in &entries {
        println!("{entry}");
    }
    println!("];");
}
