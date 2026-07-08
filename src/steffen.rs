// steffen.rs (splines library)

//! Steffen monotone cubic interpolation method.
//!
//! This is a local, shape-preserving cubic Hermite interpolant. It is usually
//! a little more conservative than PCHIP: it strongly limits node derivatives
//! around sharp changes and avoids overshoot for monotone data.

use num::Float;

use crate::alpha::cubic_coeffs_to_alpha;
use crate::binsearch::diff;
use crate::spline::Spline;

/// This function accepts x-values and y-values arrays and returns a spline interpolation container.
pub fn steffen<T: Float + std::fmt::Debug>(xx: &[T], yy: &[T]) -> Spline<T> {
    let ss = slopes_steffen(xx, yy);
    return Spline::new(xx, yy, &ss);
}

/// Estimation of tangent lines at `xx` points using Steffen's monotone method.
///
/// Endpoint slopes are the adjacent secant slopes. Interior slopes are limited
/// by Steffen's formula:
///
/// `m_i = (sign(d_{i-1}) + sign(d_i)) * min(|d_{i-1}|, |d_i|, |p_i|/2)`
///
/// where
///
/// `p_i = (d_{i-1} * h_i + d_i * h_{i-1}) / (h_{i-1} + h_i)`.
///
/// If adjacent secant slopes have opposite signs, the node slope is zero.
pub fn slopes_steffen<T: Float>(xx: &[T], yy: &[T]) -> Vec<T> {
    let n = xx.len();

    assert!(n >= 2, "steffen requires at least two points");
    assert!(yy.len() == n, "xx and yy must have equal length");

    let h: Vec<T> = diff(xx).collect();
    let delta: Vec<T> = diff(yy).zip(h.iter()).map(|(dy, dx)| dy / *dx).collect();

    if n == 2 {
        return vec![delta[0]; n];
    }

    let mut slopes = vec![T::zero(); n];
    slopes[0] = delta[0];
    slopes[n - 1] = delta[n - 2];

    for i in 1..(n - 1) {
        slopes[i] = steffen_interior_slope(h[i - 1], h[i], delta[i - 1], delta[i]);
    }

    return slopes;
}

/// Cubic coefficients for a single left Steffen interval `[x1, x2]`.
///
/// Uses the first three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn steffen_single_left<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let h12 = x2 - x1;
    let h23 = x3 - x2;
    let d12 = (y2 - y1) / h12;
    let d23 = (y3 - y2) / h23;

    let s1 = d12;
    let s2 = steffen_interior_slope(h12, h23, d12, d23);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single middle Steffen interval `[x1, x2]`.
///
/// Uses one neighbouring point on each side of the target interval. The returned
/// coefficients `[a, b, c, d]` are evaluated as `((a * dx + b) * dx + c) * dx + d`,
/// where `dx = x - x1`.
pub fn steffen_single_middle<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let h01 = x1 - x0;
    let h12 = x2 - x1;
    let h23 = x3 - x2;
    let d01 = (y1 - y0) / h01;
    let d12 = (y2 - y1) / h12;
    let d23 = (y3 - y2) / h23;

    let s1 = steffen_interior_slope(h01, h12, d01, d12);
    let s2 = steffen_interior_slope(h12, h23, d12, d23);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single right Steffen interval `[x1, x2]`.
///
/// Uses the last three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn steffen_single_right<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 4] {
    let h01 = x1 - x0;
    let h12 = x2 - x1;
    let d01 = (y1 - y0) / h01;
    let d12 = (y2 - y1) / h12;

    let s1 = steffen_interior_slope(h01, h12, d01, d12);
    let s2 = d12;

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single left Steffen interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn steffen_single_left_alpha<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = steffen_single_left(x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single middle Steffen interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn steffen_single_middle_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = steffen_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single right Steffen interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn steffen_single_right_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 2] {
    let coeffs = steffen_single_right(x0, y0, x1, y1, x2, y2);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

#[inline]
fn steffen_interior_slope<T: Float>(h1: T, h2: T, del1: T, del2: T) -> T {
    let zero = T::zero();

    if del1 * del2 <= zero {
        return zero;
    }

    let one = T::one();
    let two = one + one;
    let half = one / two;
    let p = (del1 * h2 + del2 * h1) / (h1 + h2);
    let min1 = min_abs(del1, del2);
    let min2 = min_abs(min1, half * p.abs());

    return (sign(del1) + sign(del2)) * min2;
}

#[inline]
fn sign<T: Float>(v: T) -> T {
    if v > T::zero() {
        return T::one();
    }
    if v < T::zero() {
        return -T::one();
    }
    return T::zero();
}

#[inline]
fn min_abs<T: Float>(a: T, b: T) -> T {
    if a.abs() < b.abs() {
        return a.abs();
    }
    return b.abs();
}

#[inline]
fn cubic_coeffs<T: Float>(x1: T, y1: T, x2: T, y2: T, s1: T, s2: T) -> [T; 4] {
    let dx = x2 - x1;
    let divdif = (y2 - y1) / dx;
    let dzzdx = (divdif - s1) / dx;
    let dzdxx = (s2 - divdif) / dx;
    return [(dzdxx - dzzdx) / dx, dzzdx + dzzdx - dzdxx, s1, y1];
}

#[cfg(test)]
mod tests {
    use super::{slopes_steffen, steffen_single_left, steffen_single_middle, steffen_single_right};

    fn assert_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }

    #[test]
    fn single_intervals_match_full_slope_coefficients() {
        let xx = [0.0, 1.0, 2.0, 4.0];
        let yy = [0.0, 1.0, 1.5, 3.0];
        let slopes = slopes_steffen(&xx, &yy);

        let left = steffen_single_left(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2]);
        let middle = steffen_single_middle(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);
        let right = steffen_single_right(xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);

        assert_close(left[2], slopes[0]);
        assert_close(middle[2], slopes[1]);
        assert_close(right[2], slopes[2]);
    }
}
