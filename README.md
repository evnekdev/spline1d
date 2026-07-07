# Spline1d

> A pure Rust library for 1D cubic spline interpolation, featuring Makima interpolation and tools for evaluation, root finding, and inverse lookup.

![Rust](https://img.shields.io/badge/language-Rust-orange)
![Status](https://img.shields.io/badge/status-active-green)
![License](https://img.shields.io/github/license/evnekdev/splines)

---

## 📌 Overview

Spline1d is a Rust crate for constructing and evaluating cubic spline interpolations.

### Current capabilities:
- Makima interpolation (fully implemented)
- PCHIP interpolation (work in progress)
- Piecewise polynomial evaluation (PPData)
- Fast binary search interval lookup
- Inverse interpolation via search trees
- CSV data loading

---

## ✨ Features

- Cubic spline interpolation (Makima)
- Piecewise polynomial representation
- Fast interval lookup
- Inverse interpolation
- CSV support
- Generic over floating-point types

---

## 🛠️ Tech Stack

- Rust
- num
- csv
- serde

---

## 🚀 Getting Started

git clone https://github.com/evnekdev/spline1d.git

cd splines

cargo build

cargo run

---

## 📊 Usage

```rust
use splines::makima;

fn main() {
    let x = vec![0.0, 1.0, 2.0, 3.0];
    let y = vec![0.0, 2.0, 1.0, 3.0];

    let spline = makima(&x, &y);
    let value = spline.eval(1.5);

    println!("Interpolated value: {}", value);
}
```

---

## 🧩 Core Concepts

### PPData
Stores spline coefficients and enables evaluation.

### Makima
Smooth interpolation avoiding oscillations.

### PCHIP
Shape-preserving interpolation (WIP).

### Binary Search
Efficient interval lookup.

### SearchTree
Inverse interpolation for non-monotonic curves.

### Root Solver
Analytical cubic solver.

---

## 📂 CSV Loading

use splines::load_mpp_from_csv;

---

## 🧪 Testing

cargo test

---

## ⚠️ Limitations

- Requires monotonic x values
- No NaN/Infinity handling
- PCHIP incomplete

---

## 📈 Roadmap

- Complete PCHIP
- Add 2D splines
- Improve docs
- Add benchmarks
- Visualization tools

---

## 🤝 Contributing

Fork → branch → commit → PR

---

## 🐛 Issues

https://github.com/evnekdev/splines/issues

---

## 📄 License

MIT

---

## 📬 Contact

Evgenii Nekhoroshev
https://github.com/evnekdev
