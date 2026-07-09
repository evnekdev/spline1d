# spline1d

A pure Rust library for fast 1D cubic interpolation primitives.

The crate provides both allocation-free single-interval functions and allocation-backed spline containers.

## Methods

Implemented local cubic interpolation methods:

- Akima
- Makima
- PCHIP
- Steffen
- Catmull-Rom
- Cardinal
- Fritsch-Butland

## Feature flags

```toml
[features]
default = ["std"]
std = ["alloc", "dep:csv"]
alloc = []
```

- `std` is enabled by default. It enables the full API, including CSV helpers and multi-spline/search-tree helpers.
- `alloc` enables heap-backed spline containers such as `Spline<T>` without requiring `std`.
- `--no-default-features` builds the allocation-free `no_std` API. In this mode, single-interval functions and alpha conversions remain available, but `Spline<T>`, `MultiSpline`, and CSV loading are not compiled.

## no_std usage

For allocation-free single-interval interpolation:

```toml
spline1d = { version = "0.1", default-features = false }
```

Example:

```rust
use spline1d::pchip_single_middle;

let coeffs = pchip_single_middle(
    0.0, 0.0,
    1.0, 1.0,
    2.0, 1.5,
    3.0, 2.0,
);
```

For heap-backed `Spline<T>` without `std`:

```toml
spline1d = { version = "0.1", default-features = false, features = ["alloc"] }
```

For the default desktop/server API:

```toml
spline1d = "0.1"
```

## Basic usage

```rust
use spline1d::makima;

fn main() {
    let x = vec![0.0, 1.0, 2.0, 3.0];
    let y = vec![0.0, 2.0, 1.0, 3.0];

    let spline = makima(&x, &y);
    let value = spline.interpolate(&1.5);

    println!("Interpolated value: {:?}", value);
}
```

## Current limitations

- `x` values should be monotonic.
- NaN and infinity are not supported in interval lookup.
- The `std` feature is required for CSV loading and the current `MultiSpline` / `SearchTree` helpers.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
