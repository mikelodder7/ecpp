# ECPP

[![Crate][crate-image]][crate-link]
[![Docs][docs-image]][docs-link]
![Apache-2.0/MIT licensed][license-image]
[![Downloads][downloads-image]][crate-link]
![MSRV][msrv-image]
![no_std][no-std-image]

Backend-neutral Atkin–Morain elliptic curve primality proving for Rust.

Unlike a probable-prime test, `ecpp` returns a primality certificate. The
certificate can be stored, transferred, and verified deterministically without
randomness or trust in the machine or bigint backend that produced it.

## ⚠️ Security warning

This implementation has not been independently audited.

Primality proving operates on public integers and is intentionally
variable-time. Do not use its arithmetic with secret values. A successfully
verified certificate is a mathematical proof of primality, but that does not
make an unaudited implementation immune to software defects.

**USE AT YOUR OWN RISK.**

## Installation

```toml
[dependencies]
ecpp = "0.1"
```

The default features enable `std`, operating-system randomness, and
`crypto-bigint` integration. See [Feature flags](#feature-flags) for other
configurations.

## Examples

### Prove and verify an existing prime

```rust
use ecpp::prime;
use num_bigint::BigUint;

let candidate = BigUint::from(65_537u32);
let proof = prime::prove(&candidate)?;

// Bind the proof to the expected input and verify every certificate step.
prime::verify(&candidate, &proof)?;

// A proof can also verify the number embedded in the certificate.
proof.verify()?;

# Ok::<(), ecpp::Error>(())
```

`prime::check` follows the familiar `glass_pumpkin` shape when only a boolean
result is needed:

```rust
use ecpp::prime;
use num_bigint::BigUint;

assert!(prime::check(&BigUint::from(65_537u32)));
assert!(!prime::check(&BigUint::from(561u32)));
```

### Generate a prime and its proof

`prime::new` returns both artifacts so the proof is not discarded:

```rust,no_run
use ecpp::{ProvedPrime, prime};
use num_bigint::BigUint;

let proved: ProvedPrime<BigUint> = prime::new(256)?;
proved.proof.verify_for(&proved.prime)?;

let (prime, proof) = proved.into_parts();
println!("prime: {prime}");
println!("certificate steps: {}", proof.nodes.len());

# Ok::<(), ecpp::Error>(())
```

Use `prime::from_rng` and `prime::prove_with_rng` when the application supplies
its own `rand_core::CryptoRng`.

### Use `crypto-bigint`

The `crypto-bigint` feature is enabled by default and supports both fixed-width
and boxed unsigned integers:

```rust
use crypto_bigint::U256;
use ecpp::prime;

let candidate = U256::from_u128(65_537);
let proof = prime::prove(&candidate)?;
proof.verify_for(&candidate)?;

# Ok::<(), ecpp::Error>(())
```

### Use `rug` or OpenSSL

Enable the corresponding optional feature:

```toml
[dependencies]
ecpp = { version = "0.1", features = ["rug"] }
# or: features = ["openssl"]
```

The native `rug::Integer` and `openssl::bn::BigNum` types then implement
`ecpp::Integer` and work with the same API.

### Use `unknown_order::BigNumber`

`unknown_order` can select mutually exclusive arithmetic backends. The
canonical-byte API avoids coupling `ecpp` to any particular selection:

```rust,ignore
use ecpp::prime;
use unknown_order::BigNumber;

let candidate = BigNumber::from(65_537u64);
let bytes = candidate.to_bytes();

let proof = prime::prove_be_bytes(&bytes)?;
prime::verify_be_bytes(&bytes, &proof)?;

# Ok::<(), ecpp::Error>(())
```

This path works with the `crypto`, `rust`, `gmp`, and `openssl`
`unknown_order` backends.

### Serialize a certificate

Enable `serde` to derive `Serialize` and `Deserialize` for all certificate
types. The test suite covers Postcard, CBOR via `serde_cbor_2`, JSON, TOML, and
YAML via `yaml_serde`:

```rust,ignore
use ecpp::prime;
use num_bigint::BigUint;

let candidate = BigUint::from(65_537u32);
let proof = prime::prove(&candidate)?;
let encoded = serde_json::to_vec(&proof)?;

let decoded: ecpp::PrimalityProof = serde_json::from_slice(&encoded)?;
decoded.verify_for(&candidate)?;

# Ok::<(), Box<dyn std::error::Error>>(())
```

Certificates store integers as canonical unsigned big-endian magnitudes. A
proof produced from one backend can therefore be verified against another.

## Supported backends

| Backend | Cargo feature | Integration |
| --- | --- | --- |
| `num_bigint::BigUint` / `BigInt` | Always available | Direct |
| `crypto_bigint::Uint<LIMBS>` / `BoxedUint` | `crypto-bigint` (default) | Direct |
| `rug::Integer` | `rug` | Direct |
| `openssl::bn::BigNum` | `openssl` | Direct |
| `unknown_order::BigNumber` | None required | `prove_be_bytes` / `verify_be_bytes` |
| Other bigint libraries | None required | Implement the local `ecpp::Integer` trait |

The backend trait only requires conversion to and from a canonical big-endian
magnitude. ECPP handles public values, so this conversion is deliberately
variable-time.

## Feature flags

| Feature | Default | Description |
| --- | --- | --- |
| `std` | Yes | Enables standard-library support for the arithmetic dependencies |
| `getrandom` | Yes | Enables `prime::new`, `prime::prove`, and other OS-RNG convenience functions |
| `crypto-bigint` | Yes | Implements `ecpp::Integer` for `crypto-bigint` unsigned integers |
| `rug` | No | Implements `ecpp::Integer` for `rug::Integer` |
| `openssl` | No | Implements `ecpp::Integer` for OpenSSL `BigNum` |
| `serde` | No | Enables certificate serialization and deserialization |

Without `getrandom`, use the `_with_rng` APIs. A minimal allocation-enabled
`no_std` build is available with:

```console
cargo build --no-default-features
```

## How the proof works

The prover constructs an Atkin–Morain ECPP chain. At each level it:

1. Screens the candidate for compositeness.
2. Chooses a CM discriminant and solves `4n = u² + |D|v²` with Cohen's
   modified Cornacchia algorithm.
3. Derives candidate curve orders `n + 1 ± u`.
4. Uses trial division and bounded Pollard–rho to find a sufficiently large
   probable-prime divisor `q` of a candidate order.
5. Constructs a CM elliptic curve and finds a certificate point.
6. Reduces the proof of `n` to a proof of the smaller `q`.
7. Terminates at a 64-bit prime checked by deterministic Miller–Rabin.

The CM search includes the `j = 0` and `j = 1728` supersingular cases and all
fundamental imaginary quadratic discriminants of class number one or two.

Each [`EcppStep`](https://docs.rs/ecpp/latest/ecpp/struct.EcppStep.html) records:

- the candidate `n`;
- the short Weierstrass curve `y² = x³ + ax + b (mod n)`;
- an affine point on that curve;
- a known annihilating order `m`; and
- the recursively certified prime divisor `q` of `m`.

The verifier does not trust the CM construction or the producer's
factorization work. It independently checks:

- that the curve is nonsingular modulo every divisor of `n`;
- that the certificate point lies on the curve;
- that `q` divides `m` and exceeds the elliptic Pocklington bound
  `(n^(1/4) + 1)²`;
- that multiplying by `m / q` does not annihilate the point; and
- that the following multiplication by `q` does annihilate it.

All affine divisions must be invertible modulo `n`; a nontrivial gcd causes
verification to fail.

## Search limits and errors

ECPP certificate construction is a search algorithm. This implementation uses
a bounded discriminant set, bounded factor search, and bounded point search.
Consequently, a prime may return:

```text
Error::SearchExhausted { candidate }
```

That result is inconclusive. It never means that the candidate is composite.
Prime generation handles it by trying another candidate.

`Error::Composite` means the prover found a compositeness witness during its
screening or arithmetic. `Error::InvalidProof` means a supplied certificate
failed deterministic verification.

`ProverOptions` exposes the trial-division limit, point-attempt limit, and
maximum certificate depth.

## Minimum Supported Rust Version

This crate requires **Rust 1.98** at a minimum.

The MSRV may change in a future minor release.

## Testing

The test suite covers:

- recursive construction and verification of a 128-bit certificate;
- deterministic 64-bit base cases and known pseudoprimes;
- rejection of composite inputs and tampered certificates;
- certificate serialization through Postcard, CBOR, JSON, TOML, and YAML;
- direct `num-bigint`, `crypto-bigint`, `rug`, and OpenSSL interoperability;
- canonical-byte interoperability for wrapper backends; and
- an ignored, expensive regression proving the secp256k1 field modulus.

The crate is also checked under `no_std`, documented with every independent
backend feature, and packaged with `cargo package`.

## References

The implementation and its terminology are based on the following papers and
resources:

1. Shafi Goldwasser and Joe Kilian, [*Almost All Primes Can Be Quickly
   Certified*][goldwasser-kilian-1986], Proceedings of the 18th ACM Symposium
   on Theory of Computing, 1986, pp. 316–329. This introduced elliptic-curve
   primality certificates.
2. Shafi Goldwasser and Joe Kilian, [*Primality Testing Using Elliptic
   Curves*][goldwasser-kilian-1999], Journal of the ACM 46(4), 1999,
   pp. 450–472, DOI 10.1145/320211.320213.
3. A. O. L. Atkin and François Morain, [*Elliptic Curves and Primality
   Proving*][atkin-morain], Mathematics of Computation 61(203), 1993,
   pp. 29–68, DOI 10.1090/S0025-5718-1993-1199989-X. This is the primary
   description of the practical Atkin–Morain ECPP algorithm.
4. François Morain and Jean-Louis Nicolas, [*On Cornacchia's Algorithm for
   Solving the Diophantine Equation u² + dv² = m*][cornacchia], 1990. This
   gives the norm-equation algorithm used during CM order selection.
5. Henri Cohen, [*A Course in Computational Algebraic Number
   Theory*][cohen], Graduate Texts in Mathematics 138, Springer, 1993. See the
   modified Cornacchia algorithm and computational CM background.
6. François Morain, [*Primality Proving Using Elliptic Curves: An
   Update*][morain-update], Algorithmic Number Theory, LNCS 1423, 1998,
   pp. 111–127.
7. François Morain, [*Implementing the Asymptotically Fast Version of the
   Elliptic Curve Primality Proving Algorithm*][fast-ecpp], Mathematics of
   Computation 76(257), 2007, pp. 493–505,
   DOI 10.1090/S0025-5718-06-01890-4.
8. Qi Cheng, [*Primality Proving via One Round in ECPP and One Iteration in
   AKS*][cheng], Advances in Cryptology — CRYPTO 2003, LNCS 2729,
   pp. 338–355.
9. Reinier Bröker and Peter Stevenhagen, [*Constructing Elliptic Curves with a
   Known Number of Points over a Prime Field*][broker-stevenhagen], describing
   the CM construction used in primality proving and curve generation.
10. Andrew V. Sutherland, [*Computing Hilbert Class Polynomials with the
    Chinese Remainder Theorem*][sutherland], Mathematics of Computation 80,
    2011, pp. 501–538.
11. François Morain's [ECPP home page][ecpp-home], containing implementation
    notes, historical results, certificates, and further references.

## Related projects

- [`glass_pumpkin`][glass-pumpkin] provides probable-prime testing and random
  prime generation using `num-bigint`. This crate intentionally follows its
  approachable `prime` module API while returning a proof.
- [`unknown_order`][unknown-order] provides a common API over the same family
  of multiprecision backends supported here.

## License

Licensed under either of:

- [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)
- [MIT license](http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.

[//]: # (badges)

[crate-image]: https://img.shields.io/crates/v/ecpp.svg
[crate-link]: https://crates.io/crates/ecpp
[docs-image]: https://docs.rs/ecpp/badge.svg
[docs-link]: https://docs.rs/ecpp/
[license-image]: https://img.shields.io/badge/license-Apache--2.0%2FMIT-blue.svg
[downloads-image]: https://img.shields.io/crates/d/ecpp.svg
[msrv-image]: https://img.shields.io/badge/rustc-1.98+-blue.svg
[no-std-image]: https://img.shields.io/badge/no__std-alloc-blue.svg

[atkin-morain]: https://www.ams.org/journals/mcom/1993-61-203/S0025-5718-1993-1199989-X/S0025-5718-1993-1199989-X.pdf
[broker-stevenhagen]: https://arxiv.org/abs/math/0111159
[cheng]: https://www.iacr.org/archive/crypto2003/27290337/27290337.pdf
[cohen]: https://doi.org/10.1007/978-3-662-02945-9
[cornacchia]: https://www.lix.polytechnique.fr/~morain/Articles/cornac.pdf
[ecpp-home]: https://www.lix.polytechnique.fr/~morain/Prgms/ecpp.english.html
[fast-ecpp]: https://www.lix.polytechnique.fr/~morain/Articles/fastecpp-final.pdf
[glass-pumpkin]: https://github.com/mikelodder7/glass_pumpkin
[goldwasser-kilian-1986]: https://publications.csail.mit.edu/lcs/pubs/pdf/MIT-LCS-TM-313.pdf
[goldwasser-kilian-1999]: https://doi.org/10.1145/320211.320213
[morain-update]: https://doi.org/10.1007/BFb0054855
[sutherland]: https://arxiv.org/abs/0903.2785
[unknown-order]: https://github.com/mikelodder7/unknown_order
