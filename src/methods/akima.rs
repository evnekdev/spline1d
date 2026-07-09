// akima.rs (splines library)

//! Original Akima interpolation method.
//!
//! This is the unmodified Akima local cubic Hermite interpolant. Node slopes
//! are calculated from four neighboring secant slopes using Akima's weighted
//! formula. Missing endpoint secant slopes are obtained by linear extrapolation,
//! matching the endpoint treatment used by the existing Makima implementation.

use num_traits::Float;

#[cfg(feature = "alloc")]
use alloc::vec;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

use crate::alpha::cubic_coeffs_to_alpha;
use crate::binsearch::diff;
#[cfg(feature = "alloc")]
use crate::spline::Spline;

/// This function accepts x-values and y-values arrays and returns a spline interpolation container.
#[cfg(feature = "alloc")]
pub fn akima<T: Float + core::fmt::Debug>(xx: &[T], yy: &[T]) -> Spline<T> {
    let ss = slopes_akima(xx, yy);
    return Spline::new(xx, yy, &ss);
}

/// Estimation of tangent lines at `xx` points using Akima's original method.
///
/// For interval secant slopes `delta_i = (y_{i+1} - y_i)/(x_{i+1} - x_i)`,
/// the slope at node `i` is
///
/// `d_i = (w1 * delta_{i-1} + w2 * delta_i) / (w1 + w2)`,
///
/// where `w1 = |delta_{i+1} - delta_i|` and
/// `w2 = |delta_{i-1} - delta_{i-2}|`.
///
/// If both weights are zero, the conventional Akima fallback
/// `(delta_{i-1} + delta_i)/2` is used.
#[cfg(feature = "alloc")]
pub fn slopes_akima<T: Float>(xx: &[T], yy: &[T]) -> Vec<T> {
    let n = xx.len();

    assert!(n >= 2, "akima requires at least two points");
    assert!(yy.len() == n, "xx and yy must have equal length");

    let hh: Vec<T> = diff(xx).collect();
    let delta: Vec<T> = diff(yy).zip(hh.iter()).map(|(dy, dx)| dy / *dx).collect();

    // Special case n = 2, use linear interpolation.
    if n == 2 {
        return vec![delta[0]; n];
    }

    // Akima requires two additional secant slopes on each side. Following the
    // common Akima/Makima endpoint convention, extrapolate those slopes linearly.
    let delta_m1 = (delta[0] + delta[0]) - delta[1];
    let delta_m2 = (delta_m1 + delta_m1) - delta[0];
    let delta_n = (delta[n - 2] + delta[n - 2]) - delta[n - 3];
    let delta_n1 = (delta_n + delta_n) - delta[n - 2];

    let mut ext = Vec::with_capacity(n + 3);
    ext.push(delta_m2);
    ext.push(delta_m1);
    ext.extend(delta.iter().copied());
    ext.push(delta_n);
    ext.push(delta_n1);

    let mut slopes = Vec::with_capacity(n);
    for i in 0..n {
        slopes.push(akima_slope(ext[i], ext[i + 1], ext[i + 2], ext[i + 3]));
    }

    return slopes;
}

/// Cubic coefficients for a single left Akima interval `[x1, x2]`.
///
/// Uses only the first three data points and Akima's linear extrapolation of
/// missing secant slopes on the left. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn akima_single_left<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let d12 = secant(x1, y1, x2, y2);
    let d23 = secant(x2, y2, x3, y3);

    let dm1 = (d12 + d12) - d23;
    let dm2 = (dm1 + dm1) - d12;
    let d34 = (d23 + d23) - d12;

    let s1 = akima_slope(dm2, dm1, d12, d23);
    let s2 = akima_slope(dm1, d12, d23, d34);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single middle Akima interval `[x1, x2]`.
///
/// Uses one neighbouring point on each side of the target interval. The returned
/// coefficients `[a, b, c, d]` are evaluated as `((a * dx + b) * dx + c) * dx + d`,
/// where `dx = x - x1`.
pub fn akima_single_middle<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let d01 = secant(x0, y0, x1, y1);
    let d12 = secant(x1, y1, x2, y2);
    let d23 = secant(x2, y2, x3, y3);

    let dm1 = (d01 + d01) - d12;
    let d34 = (d23 + d23) - d12;

    let s1 = akima_slope(dm1, d01, d12, d23);
    let s2 = akima_slope(d01, d12, d23, d34);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single right Akima interval `[x1, x2]`.
///
/// Uses only the last three data points and Akima's linear extrapolation of
/// missing secant slopes on the right. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn akima_single_right<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 4] {
    let d01 = secant(x0, y0, x1, y1);
    let d12 = secant(x1, y1, x2, y2);

    let dm1 = (d01 + d01) - d12;
    let d23 = (d12 + d12) - d01;
    let d34 = (d23 + d23) - d12;

    let s1 = akima_slope(dm1, d01, d12, d23);
    let s2 = akima_slope(d01, d12, d23, d34);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single left Akima interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn akima_single_left_alpha<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = akima_single_left(x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single middle Akima interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn akima_single_middle_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = akima_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single right Akima interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn akima_single_right_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 2] {
    let coeffs = akima_single_right(x0, y0, x1, y1, x2, y2);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

#[inline]
fn secant<T: Float>(xa: T, ya: T, xb: T, yb: T) -> T {
    return (yb - ya) / (xb - xa);
}

#[inline]
fn akima_slope<T: Float>(dm2: T, dm1: T, d0: T, d1: T) -> T {
    let w1 = (d1 - d0).abs();
    let w2 = (dm1 - dm2).abs();
    let w = w1 + w2;

    if w == T::zero() {
        let half = T::one() / (T::one() + T::one());
        return (dm1 + d0) * half;
    }

    return (w1 * dm1 + w2 * d0) / w;
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
    use super::{akima_single_left, akima_single_middle, akima_single_right, slopes_akima};

    fn assert_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }

    #[test]
    fn single_intervals_match_full_slope_coefficients() {
        let xx = [0.0, 1.0, 2.0, 4.0];
        let yy = [0.0, 1.0, 1.5, 3.0];
        let slopes = slopes_akima(&xx, &yy);

        let left = akima_single_left(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2]);
        let middle = akima_single_middle(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);
        let right = akima_single_right(xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);

        assert_close(left[2], slopes[0]);
        assert_close(middle[2], slopes[1]);
        assert_close(right[2], slopes[2]);
    }
}
