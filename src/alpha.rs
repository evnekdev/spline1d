// alpha.rs (splines library)

//! Conversion utilities for the normalized endpoint form.
//!
//! `y(t) = y1 * (1 - t) + y2 * t + (1 - t) * t * (alpha0 + alpha1 * t)`,
//! where `t = (x - x1) / (x2 - x1)`.

use num_traits::Float;

/// Convert direct cubic coefficients `[a, b, c, d]` on `dx = x - x1` to
/// normalized interval coefficients `[alpha0, alpha1]`.
///
/// The direct cubic is evaluated as:
///
/// `y = ((a * dx + b) * dx + c) * dx + d`
///
/// The normalized form is:
///
/// `y = y1 * (1 - t) + y2 * t + (1 - t) * t * (alpha0 + alpha1 * t)`
///
/// with `t = dx / h` and `h = x2 - x1`.
///
/// Note that the interval width `h` is required. The four direct coefficients
/// alone are not enough unless the interval width is known or assumed to be 1.
#[inline]
pub fn cubic_coeffs_to_alpha<T: Float>(coeffs: [T; 4], h: T) -> [T; 2] {
    let a = coeffs[0];
    let b = coeffs[1];
    let h2 = h * h;
    let h3 = h2 * h;

    let alpha1 = -a * h3;
    let alpha0 = alpha1 - b * h2;

    return [alpha0, alpha1];
}

/// Convert direct cubic coefficients to `[alpha0, alpha1]` for a unit-width
/// interval. This is equivalent to `cubic_coeffs_to_alpha(coeffs, T::one())`.
#[inline]
pub fn cubic_coeffs_to_alpha_unit<T: Float>(coeffs: [T; 4]) -> [T; 2] {
    return cubic_coeffs_to_alpha(coeffs, T::one());
}

/// Convert normalized interval coefficients `[alpha0, alpha1]` back to direct
/// cubic coefficients `[a, b, c, d]` on `dx = x - x1`.
///
/// The returned coefficients are evaluated as:
///
/// `y = ((a * dx + b) * dx + c) * dx + d`
///
/// where `dx = x - x1`.
#[inline]
pub fn alpha_to_cubic_coeffs<T: Float>(x1: T, y1: T, x2: T, y2: T, alpha0: T, alpha1: T) -> [T; 4] {
    let h = x2 - x1;
    let dy = y2 - y1;
    let h2 = h * h;
    let h3 = h2 * h;

    let a = -alpha1 / h3;
    let b = (alpha1 - alpha0) / h2;
    let c = (dy + alpha0) / h;
    let d = y1;

    return [a, b, c, d];
}

/// Convert normalized interval coefficients `[alpha0, alpha1]` to standard
/// global cubic coefficients `[A, B, C, D]` in `x`:
///
/// `y = ((A * x + B) * x + C) * x + D`.
///
/// Prefer `alpha_to_cubic_coeffs` for spline evaluation, because the local
/// `dx = x - x1` form is usually numerically better conditioned.
#[inline]
pub fn alpha_to_standard_cubic_coeffs<T: Float>(x1: T, y1: T, x2: T, y2: T, alpha0: T, alpha1: T) -> [T; 4] {
    let local = alpha_to_cubic_coeffs(x1, y1, x2, y2, alpha0, alpha1);
    let a = local[0];
    let b = local[1];
    let c = local[2];
    let d = local[3];

    let two = T::one() + T::one();
    let three = two + T::one();
    let x1_2 = x1 * x1;
    let x1_3 = x1_2 * x1;

    let aa = a;
    let bb = b - three * a * x1;
    let cc = c - two * b * x1 + three * a * x1_2;
    let dd = d - c * x1 + b * x1_2 - a * x1_3;

    return [aa, bb, cc, dd];
}

#[cfg(test)]
mod tests {
    use super::{alpha_to_cubic_coeffs, cubic_coeffs_to_alpha};

    fn eval_local(coeffs: [f64; 4], dx: f64) -> f64 {
        return ((coeffs[0] * dx + coeffs[1]) * dx + coeffs[2]) * dx + coeffs[3];
    }

    fn eval_alpha(x1: f64, y1: f64, x2: f64, y2: f64, alpha0: f64, alpha1: f64, x: f64) -> f64 {
        let t = (x - x1) / (x2 - x1);
        return y1 * (1.0 - t) + y2 * t + (1.0 - t) * t * (alpha0 + alpha1 * t);
    }

    #[test]
    fn alpha_roundtrip_matches_local_cubic() {
        let x1 = 2.0;
        let x2 = 5.0;
        let h = x2 - x1;
        let coeffs = [0.25, -1.5, 3.0, 7.0];
        let y1 = coeffs[3];
        let y2 = eval_local(coeffs, h);
        let [alpha0, alpha1] = cubic_coeffs_to_alpha(coeffs, h);
        let roundtrip = alpha_to_cubic_coeffs(x1, y1, x2, y2, alpha0, alpha1);

        for k in 0..4 {
            assert!((coeffs[k] - roundtrip[k]).abs() < 1e-12);
        }

        for x in [2.0, 2.75, 3.5, 4.25, 5.0] {
            let y_cubic = eval_local(coeffs, x - x1);
            let y_alpha = eval_alpha(x1, y1, x2, y2, alpha0, alpha1, x);
            assert!((y_cubic - y_alpha).abs() < 1e-12);
        }
    }
}
