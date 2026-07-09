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

## Inverse lookup and multi-value interpolation

`spline1d` supports inverse lookup through `SearchTree` and `SearchNode`.

This is useful when the independent and dependent variables need to be swapped, for example when finding all values of `x` corresponding to a target value of `y`.

Unlike a simple binary search, inverse lookup cannot assume that the spline is globally monotone. A single value of `y` may correspond to several different `x` values. `SearchTree` handles this by splitting the spline into monotone regions and searching only the regions whose value ranges contain the requested target.

Conceptually, this allows queries such as:

```text
x -> y
```

and also:

```text
y -> all matching x values
```

The tree is built from a `MultiSpline`:

```rust
use spline1d::*;

let tree = SearchTree::new(&multispline);
```

Then interpolation can be performed between any registered variables:

```rust
let values = tree.interpolate(&key_x, &key_y, &value);
```

Depending on the selected keys, this can perform:

* ordinary interpolation from the principal coordinate to a dependent variable;
* inverse interpolation from a dependent variable back to the principal coordinate;
* cross-variable interpolation through the principal coordinate;
* identity lookup when both keys refer to the principal coordinate.

For non-monotone data, the returned `Vec<T>` may contain multiple values. This means `spline1d` can find all valid solutions instead of returning only the first match or requiring the user to scan intervals manually.

Internally, each `SearchNode` stores an interval of spline indices and value ranges for every variable. During lookup, branches whose min/max range cannot contain the target value are skipped. Candidate monotone regions are then searched and evaluated locally.


## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
