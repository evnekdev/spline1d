//! Finding the cubic root (y = 0 intersection) analytically.
#[allow(unused_imports)]

use num_traits::Float;

fn powf_(x: f64, pow: f64) -> f64 {
    let absroot = x.abs().powf(pow);
    if x < 0.0 {
        return -absroot;
    }
    return absroot;
}

/// Calculate a cubic root to find x intersection of a cubic spline
pub fn calculate_root(coeffs: &[f64]) -> f64 {
    let a = coeffs[0];
    let b = coeffs[1];
    let c = coeffs[2];
    let d = coeffs[3];
    let q = (3.0 * a * c - b * b) / (9.0 * a * a);
    let r = (9.0 * a * b * c - 27.0 * a * a * d - 2.0 * b * b * b) / (54.0 * a * a * a);
    let sqrtvalue = (q * q * q + r * r).sqrt();
    let s = powf_(r + sqrtvalue, 1.0 / 3.0);
    let t = powf_(r - sqrtvalue, 1.0 / 3.0);
    let x1 = s + t - b / (3.0 * a);
    return x1;
}
