# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
