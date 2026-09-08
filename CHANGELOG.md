# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-09-08

### Added

- The CM discriminant table grew from class numbers one and two (25
  entries) to every fundamental discriminant with class number up to eight
  (333 entries), generated with verified high-precision arithmetic and
  stored as exact byte-encoded coefficients. A randomized recursive
  Cantor–Zassenhaus splitter finds roots of the higher-degree Hilbert
  class polynomials in both the generic engine and the allocation-free
  `fixed` module.
- Backtracking proof search: a failed level no longer abandons the proof.
  Each level tries every available curve order, re-sweeps with fresh
  randomness, and backtracks to alternate branches, bounded by the new
  `ProverOptions::search_budget`. `Error::SearchExhausted` now reports the
  original candidate rather than an internal intermediate.
- Together, the enlarged table, backtracking, and ECM move the practical
  proving envelope from roughly 512 bits to 2048 bits: a 1024-bit prime
  that previously exhausted the search in minutes now proves with default
  options in tens of minutes, and a 2048-bit prime proves overnight
  (one seeded search, run twice: 14.5 and 15.2 hours for the same 37-step
  certificate, with under a fifth of the default search budget consumed).

- Lenstra ECM (Montgomery x-only stage one with Suyama parameterization) as
  an escalating fallback when Pollard rho cannot split a curve-order
  cofactor of at least 384 bits, in both the generic engine and the
  allocation-free `fixed` module. This extends the reachable factor range
  well beyond rho.
- README guidance on when to use ECPP over probabilistic testing, with
  measured proving, verification, and screening costs across input sizes.

### Changed

- **Breaking:** both `ProverOptions` types gained three fields —
  `max_class_number` bounding which CM discriminants are searched,
  `ecm_rounds` controlling ECM escalation, and `search_budget` bounding
  total backtracking work — so struct-literal constructions must add them
  (or use `..Default::default()`). The heap prover defaults to the full
  table (class number eight), one ECM round, a budget of 8192, and a
  `max_depth` of 128 so deep descents are not cut short. The
  allocation-free prover keeps conservative defaults (class number two, no
  ECM, budget 1024) so constrained `no_std` targets never pay for the
  heavier search unless they opt in.

## [0.1.1] - 2026-09-03

### Added

- `ArithmeticBackend::sqrt` provided method (integer square root, rounding
  down) with a Newton default and native implementations for `num-bigint`,
  `rug`, and both `crypto-bigint` backends. The engine's hand-rolled
  `integer_sqrt` was removed in its favor.

### Changed

- The heap-backed `crypto-bigint` backend (`CryptoBigint`) now delegates
  `modular_pow`, `gcd`, and `modular_inverse` to native `BoxedUint`
  implementations (Montgomery exponentiation, safegcd, and `invert_mod`)
  instead of the generic operator-based algorithms. The 128-bit
  `crypto-bigint` proving test runs about 4.6 times faster.

### Fixed

- `CryptoBigint::modular_inverse` returned `Error::Composite` for a modulus
  of one; it now returns zero, matching every other backend.
- The fixed-width `crypto-bigint` backend no longer delegates `jacobi` to
  `Uint::jacobi_symbol_vartime`, which returns an incorrect sign for some
  inputs of four or more limbs in every released `crypto-bigint` through
  0.7.5 (RustCrypto/crypto-bigint#1295, fixed upstream but unreleased). The
  generic binary algorithm is used instead.
- The allocation-free prover in `fixed` computes Jacobi symbols with an
  in-module binary algorithm instead of the affected
  `Uint::jacobi_symbol_vartime`. The wrong sign could make the prover miss
  usable curves or square roots; certificate verification was never
  affected.

## [0.1.0] - 2026-09-02

### Added

- Initial release: backend-neutral Atkin–Morain ECPP proving and
  verification over `num-bigint`, `crypto-bigint` (fixed-width and
  heap-backed), `rug`, and OpenSSL, plus the allocation-free `fixed` module
  for `no_std` targets without an allocator.
