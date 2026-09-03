# ECPP

[![Crate][crate-image]][crate-link]
[![Docs][docs-image]][docs-link]
![Apache-2.0/MIT licensed][license-image]
[![Downloads][downloads-image]][crate-link]
![MSRV][msrv-image]
![no_std][no-std-image]

Backend-neutral Atkin–Morain Elliptic Curve Primality Proving.

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

The default features enable allocation, `std`, operating-system randomness,
and the `num-bigint` and `crypto-bigint` backends. See
[Feature flags](#feature-flags) for other configurations.

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
use ecpp::{PrimalityProof, ProvedPrime, prime};
use num_bigint::BigUint;

let proved: ProvedPrime<BigUint, PrimalityProof> = prime::new(256)?;
proved.proof.verify_for(&proved.prime)?;

let (prime, proof) = proved.into_parts();
println!("prime: {prime}");
println!("certificate steps: {}", proof.nodes.len());

# Ok::<(), ecpp::Error>(())
```

Use `prime::from_rng` and `prime::prove_with_rng` when the application supplies
its own `rand_core::CryptoRng`.

### Use `crypto-bigint` without allocation

The `crypto-bigint` feature is enabled by default. Its `fixed` module performs
the complete proof with `Uint<LIMBS>` and a caller-selected stack capacity:

```rust
use crypto_bigint::U256;
use ecpp::fixed::{PrimalityProof, prove_with_rng};
use rand::{SeedableRng, rngs::StdRng};

let candidate = U256::from_u128(65_537);
let mut rng = StdRng::seed_from_u64(7);
let proof: PrimalityProof<{ U256::BYTES }, 64> =
    prove_with_rng(&candidate, &mut rng)?;
proof.verify_for(&candidate)?;

# Ok::<(), ecpp::fixed::Error>(())
```

### Use `rug` or OpenSSL

Enable the corresponding optional feature:

```toml
[dependencies]
ecpp = { version = "*", default-features = false, features = ["getrandom", "rug"] }
# or: features = ["getrandom", "openssl"]
```

The native `rug::Integer` and `openssl::bn::BigNum` types then implement
`ecpp::Integer` and work with the same API. Each feature supplies its own
arithmetic engine; neither enables `num-bigint`.

```rust
use ecpp::prime;
use rug::Integer;

let candidate = Integer::from(65_537);
let proof = prime::prove(&candidate)?;
proof.verify_for(&candidate)?;

# Ok::<(), ecpp::Error>(())
```

When several engines are enabled together, `prime::prove_with_backend` and
`PrimalityProof::verify_with` let the caller choose explicitly. This also makes
cross-backend verification visible at the call site:

```rust,ignore
let proof = prime::prove_with_backend::<ecpp::arithmetic::Rug, _, _>(
    &candidate,
    &mut rng,
    ecpp::ProverOptions::default(),
)?;
proof.verify_with::<ecpp::arithmetic::OpenSsl>()?;
```

### Serialize a certificate

Enable `serde` to implement `Serialize` and `Deserialize` for all certificate
types. The test suite covers `postcard`, CBOR, JSON, TOML, and YAML:

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

Certificates store integers as canonical unsigned big-endian magnitudes.
Human-readable formats serialize those magnitudes as lowercase hexadecimal
strings, while binary formats use byte strings. A proof produced from one
backend can therefore be verified against another.

## Supported backends

| Backend | Cargo feature | Integration |
| --- | --- | --- |
| `num_bigint::BigUint` / `BigInt` | `num-bigint` (default) | Direct |
| `crypto_bigint::Uint<LIMBS>` | `crypto-bigint` (default) | Allocation-free fixed proof; generic `CryptoUint<LIMBS>` with `alloc` |
| `crypto_bigint::BoxedUint` | `crypto-bigint`, `alloc` | Generic `CryptoBigint` engine |
| `rug::Integer` | `rug` | Direct |
| `openssl::bn::BigNum` | `openssl` | Direct |
| `unknown_order::BigNumber` | Any heap engine | `prove_be_bytes` / `verify_be_bytes` |
| Other bigint libraries | Any heap engine | Implement `ecpp::Integer`, or use canonical bytes |

`Integer` only converts caller values to and from a canonical big-endian
magnitude. `ArithmeticBackend` drives the generic ECPP engine. The default
engine is `num-bigint` when enabled, followed by boxed `crypto-bigint`, `rug`,
and OpenSSL. The `*_with_backend` APIs select one explicitly. ECPP handles
public values, so these operations are deliberately variable-time.

## Feature flags

| Feature | Default | Description |
| --- | --- | --- |
| `alloc` | Yes | Enables neutral heap proof types and generic APIs without selecting an engine |
| `std` | Yes | Extends enabled backends with standard-library support |
| `num-bigint` | Yes | Enables the `num-bigint` arithmetic engine |
| `getrandom` | Yes | Enables `prime::new`, `prime::prove`, and other OS-RNG convenience functions |
| `crypto-bigint` | Yes | Enables allocation-free fixed proving and, with `alloc`, boxed and fixed-width generic engines |
| `rug` | No | Implements `ecpp::Integer` for `rug::Integer` |
| `openssl` | No | Implements `ecpp::Integer` for OpenSSL `BigNum` |
| `serde` | No | Enables certificate serialization and deserialization |

Without `getrandom`, use the `_with_rng` APIs. A `num-bigint`-backed `no_std`
build is available with:

```console
cargo build --no-default-features --features num-bigint
```

Fixed-width `crypto-bigint` proving needs no allocator and does not activate
`num-bigint`:

```console
cargo build --no-default-features --features crypto-bigint
```

```rust
use crypto_bigint::U256;
use ecpp::fixed::{PrimalityProof, from_rng, prove_with_rng};
use ecpp::ProvedPrime;
use rand::{SeedableRng, rngs::StdRng};

let candidate = U256::from(65_537u32);
let mut rng = StdRng::seed_from_u64(7);
let proof: PrimalityProof<{ U256::BYTES }, 64> =
    prove_with_rng(&candidate, &mut rng)?;
proof.verify_for(&candidate)?;

let generated: ProvedPrime<U256, PrimalityProof<{ U256::BYTES }, 64>> =
    from_rng(128, &mut rng)?;
generated.proof.verify_for(&generated.prime)?;

# Ok::<(), Box<dyn std::error::Error>>(())
```

When `alloc` is also enabled, `proof.to_alloc()` converts a fixed certificate
to the neutral heap format. That certificate can then be verified using any
generic engine. `arithmetic::CryptoBigint` is backed by `BoxedUint` and grows
as needed. `arithmetic::CryptoUint<LIMBS>` uses `Uint<LIMBS>` as an explicitly
bounded working representation:

```rust
use crypto_bigint::{U256, U512};
use ecpp::{ProverOptions, arithmetic::CryptoUint, prime};
use rand::{SeedableRng, rngs::StdRng};

let candidate = U256::from(65_537u32);
let mut rng = StdRng::seed_from_u64(7);
let proof = prime::prove_with_backend::<CryptoUint<{ U512::LIMBS }>, _, _>(
    &candidate,
    &mut rng,
    ProverOptions::default(),
)?;
proof.verify_for_with::<CryptoUint<{ U512::LIMBS }>, _>(&candidate)?;

# Ok::<(), ecpp::Error>(())
```

The working width should normally be at least twice the candidate width so
unreduced products fit. Insufficient capacity returns `Error::Arithmetic`
instead of wrapping.

## How the proof works

The prover constructs an Atkin–Morain ECPP chain. At each level it:

1. Screens the candidate for compositeness.
2. Chooses a CM discriminant and solves `4n = u² + |D|v²` with Cohen's
   modified Cornacchia algorithm.
3. Derives candidate curve orders `n + 1 ± u`.
4. Uses trial division and bounded Pollard rho to find a sufficiently large
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
- the cofactor `m/q`; and
- the recursively certified prime divisor `q`.

The verifier does not trust the CM construction or the producer's
factorization work. It independently checks:

- that the curve is nonsingular modulo every divisor of `n`;
- that the certificate point lies on the curve;
- that `q` exceeds the elliptic Pocklington bound
  `(n^(1/4) + 1)²`;
- that multiplying by the recorded cofactor `m/q` does not annihilate the
  point; and
- that the following multiplication by `q` does annihilate it.

All affine divisions must be invertible modulo `n`; a nontrivial gcd causes
verification to fail.

## Search limits and errors

ECPP certificate construction is a search algorithm. This implementation uses
a bounded discriminant set, bounded factor search, and bounded point search.
Consequently, proving a prime may return:

```text
Error::SearchExhausted { candidate }
```

That result is inconclusive. It never means that the candidate is composite.
Prime generation handles it by trying another candidate.

`Error::Composite` means the prover found a compositeness witness during its
screening or arithmetic. `Error::InvalidProof` means a supplied certificate
failed deterministic verification.

`ProverOptions` exposes the trial-division limit, point-search limit, and
maximum certificate depth.

## Minimum supported Rust version

This crate requires **Rust 1.98** at a minimum.

The MSRV may change in a future minor release.

## Testing

The test suite covers:

- recursive construction and verification of a 128-bit certificate;
- deterministic 64-bit base cases and known pseudoprimes;
- rejection of composite inputs and tampered certificates;
- certificate serialization through `postcard`, CBOR, JSON, TOML, and YAML;
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

- [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0)
- [MIT license](https://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this work by you, as defined in the Apache-2.0
license, shall be dual-licensed as above, without any additional terms or
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
