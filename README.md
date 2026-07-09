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

## Inverse lookup

In addition to conventional interpolation (`x → y`), `spline1d` provides efficient **inverse lookup** (`y → x`).

Given a target value `y`, the library can locate **every** corresponding `x` value satisfying

```text
Spline(x) = y
```

This capability is particularly useful for calibration curves, lookup tables, engineering property databases, phase diagrams, signal processing, and other applications where a function must be queried in both directions.

Inverse lookup is implemented using the `SearchTree` and `SearchNode` structures. Rather than scanning every spline interval, a search tree organizes intervals according to their value ranges, allowing only potentially matching intervals to be examined before solving the corresponding cubic equations.

```rust
use spline1d::prelude::*;

let xs = vec![0.0, 1.0, 2.0, 3.0];
let ys = vec![0.0, 2.0, 1.0, 4.0];

let spline = Spline::makima(&xs, &ys);

// Construct a search tree for efficient inverse lookup.
let tree = SearchTree::new(&spline);

// Find every x such that spline(x) = 1.5.
let xs = tree.inverse(1.5);
```

Unlike many interpolation libraries that assume the interpolated function is monotone and return at most one solution, `spline1d` supports **non-monotone** splines and returns **all** valid solutions. This makes it suitable for functions with multiple local extrema, where the same output value may occur at several distinct locations.

The search tree is reusable, making repeated inverse queries efficient even for large splines.


## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
