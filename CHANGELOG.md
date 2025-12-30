# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2025-12-30

### Changed

- **Rebranding to Luminescent**: Renamed package from `curvelib` to `lumen-curve`
  - Now part of the [luminescent](https://github.com/Tranduy1dol/luminescent) project
  - Updated all import paths: `curvelib::*` → `lumen_curve::*`
  - Dependency updated from `mathlib` to `lumen-math`
  - Repository URL changed to `https://github.com/Tranduy1dol/lumen-curve`

### Migration Guide

```rust
// Old (curvelib)
use curvelib::instances::tiny_jubjub::TinyJubjubConfig;
use curvelib::protocol::signing::SigningEngine;

// New (lumen-curve)
use lumen_curve::instances::tiny_jubjub::TinyJubjubConfig;
use lumen_curve::protocol::signing::SigningEngine;
```

## [0.1.0] - 2025-12-15

### Added

- **Elliptic Curve Models**:
  - **Short Weierstrass**: Generic implementation for curves like secp256k1
  - **Twisted Edwards**: Support for Edwards curves like Jubjub
  - **Sextic Twist**: Implementation of twist curves for pairing-friendly constructions

- **Pairing-Friendly Curves**:
  - Implementation of bilinear pairings (Miller loop)
  - Support for BN and BLS curve families

- **Cryptographic Primitives**:
  - **KeyEngine**: Secure keypair generation and serialization
  - **SigningEngine**: ECDSA-like signing and verification
  - **Pedersen Commitments**: (in progress)

- **Curve Instances**:
  - **Tiny Jubjub**: A toy example for educational purposes
  - **BLS6_6**: A small embedding degree curve for testing pairings

- **Core Traits**:
  - `CurveConfig`: Define curve parameters
  - `GroupElement`: Point operations
  - `Pairing`: Bilinear pairing interface
