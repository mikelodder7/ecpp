# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Twelve class-number-three CM discriminants (`-23` through `-547`) with
  exact cubic Hilbert class polynomials, plus a randomized Cantor–Zassenhaus
  cubic root-finder in both the generic engine and the allocation-free
  `fixed` module. The enlarged search proves candidates that previously
  returned `Error::SearchExhausted`; a 512-bit prime that exhausted the old
  table now proves with default `ProverOptions`.
- README guidance on when to use ECPP over probabilistic testing, with
  measured proving, verification, and screening costs across input sizes.

- Lenstra ECM (Montgomery x-only stage one with Suyama parameterization) as
  an escalating fallback when Pollard rho cannot split a curve-order
  cofactor of at least 384 bits, in both the generic engine and the
  allocation-free `fixed` module. This extends the reachable factor range
  well beyond rho and is the first step toward proving 1024-bit candidates.

### Changed

- Both `ProverOptions` types gained an `ecm_rounds` field controlling ECM
  escalation (heap prover defaults to one round; the allocation-free prover
  defaults to zero so constrained targets never pay for it), and a
  `max_class_number` field bounding
  which CM discriminants are searched. The heap prover defaults to three
  (the full table); the allocation-free prover defaults to two, so
  constrained `no_std` targets skip cubic polynomial splitting unless they
  opt in. Struct-literal constructions of `ProverOptions` must add the new
  field.

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
