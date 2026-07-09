// fritschbutland.rs (splines library)

//! Fritsch-Butland monotone cubic interpolation method.
//!
//! This module implements the symmetric local nonlinear averaging function
//! described by Fritsch and Butland before their h-dependent Brodlie/PCHIP
//! modification.  It is a monotone cubic Hermite method: if two adjacent secant
//! slopes have opposite signs or either is zero, the node derivative is zero;
//! otherwise a limited harmonic-type mean is used.

use num::Float;

use crate::alpha::cubic_coeffs_to_alpha;
use crate::binsearch::diff;
use crate::spline::Spline;

/// This function accepts x-values and y-values arrays and returns a spline interpolation container.
pub fn fritschbutland<T: Float + std::fmt::Debug>(xx: &[T], yy: &[T]) -> Spline<T> {
    let ss = slopes_fritschbutland(xx, yy);
    return Spline::new(xx, yy, &ss);
}

/// Estimation of tangent lines at `xx` points using the Fritsch-Butland limiter.
///
/// Interior slopes are calculated from adjacent secant slopes `s1` and `s2` as
///
/// `0`, if `s1 * s2 <= 0`, otherwise
///
/// `sign(s1) * 3 * abs(s1) * abs(s2) / (max(abs(s1), abs(s2)) + 2 * min(abs(s1), abs(s2)))`.
///
/// Endpoint slopes use the same shape-preserving one-sided three-point formula
/// used by PCHIP/Netlib PCHIM.
pub fn slopes_fritschbutland<T: Float>(xx: &[T], yy: &[T]) -> Vec<T> {
    let n = xx.len();

    assert!(n >= 2, "fritschbutland requires at least two points");
    assert!(yy.len() == n, "xx and yy must have equal length");

    let h: Vec<T> = diff(xx).collect();
    let delta: Vec<T> = diff(yy).zip(h.iter()).map(|(dy, dx)| dy / *dx).collect();

    if n == 2 {
        return vec![delta[0]; n];
    }

    let mut slopes = vec![T::zero(); n];
    slopes[0] = endpoint_slope(h[0], h[1], delta[0], delta[1]);
    slopes[n - 1] = endpoint_slope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);

    for i in 1..(n - 1) {
        slopes[i] = fritschbutland_interior_slope(delta[i - 1], delta[i]);
    }

    return slopes;
}

/// Cubic coefficients for a single left Fritsch-Butland interval `[x1, x2]`.
///
/// Uses the first three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn fritschbutland_single_left<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let h12 = x2 - x1;
    let h23 = x3 - x2;
    let d12 = (y2 - y1) / h12;
    let d23 = (y3 - y2) / h23;

    let s1 = endpoint_slope(h12, h23, d12, d23);
    let s2 = fritschbutland_interior_slope(d12, d23);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single middle Fritsch-Butland interval `[x1, x2]`.
///
/// Uses one neighbouring point on each side of the target interval. The returned
/// coefficients `[a, b, c, d]` are evaluated as `((a * dx + b) * dx + c) * dx + d`,
/// where `dx = x - x1`.
pub fn fritschbutland_single_middle<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let d01 = secant(x0, y0, x1, y1);
    let d12 = secant(x1, y1, x2, y2);
    let d23 = secant(x2, y2, x3, y3);

    let s1 = fritschbutland_interior_slope(d01, d12);
    let s2 = fritschbutland_interior_slope(d12, d23);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single right Fritsch-Butland interval `[x1, x2]`.
///
/// Uses the last three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn fritschbutland_single_right<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 4] {
    let h01 = x1 - x0;
    let h12 = x2 - x1;
    let d01 = (y1 - y0) / h01;
    let d12 = (y2 - y1) / h12;

    let s1 = fritschbutland_interior_slope(d01, d12);
    let s2 = endpoint_slope(h12, h01, d12, d01);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single left Fritsch-Butland interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn fritschbutland_single_left_alpha<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = fritschbutland_single_left(x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single middle Fritsch-Butland interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn fritschbutland_single_middle_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = fritschbutland_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single right Fritsch-Butland interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn fritschbutland_single_right_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 2] {
    let coeffs = fritschbutland_single_right(x0, y0, x1, y1, x2, y2);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

#[inline]
fn secant<T: Float>(xa: T, ya: T, xb: T, yb: T) -> T {
    return (yb - ya) / (xb - xa);
}

#[inline]
fn fritschbutland_interior_slope<T: Float>(del1: T, del2: T) -> T {
    let zero = T::zero();

    if del1 * del2 <= zero {
        return zero;
    }

    let one = T::one();
    let two = one + one;
    let three = two + one;
    let abs1 = del1.abs();
    let abs2 = del2.abs();
    let dmax = if abs1 > abs2 { abs1 } else { abs2 };
    let dmin = if abs1 < abs2 { abs1 } else { abs2 };

    return sign(del1) * three * abs1 * abs2 / (dmax + two * dmin);
}

#[inline]
fn endpoint_slope<T: Float>(h1: T, h2: T, del1: T, del2: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let two = one + one;
    let three = two + one;

    let mut d = ((two * h1 + h2) * del1 - h1 * del2) / (h1 + h2);

    if d * del1 <= zero {
        d = zero;
    } else if del1 * del2 < zero && d.abs() > (three * del1).abs() {
        d = three * del1;
    }

    return d;
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
fn cubic_coeffs<T: Float>(x1: T, y1: T, x2: T, y2: T, s1: T, s2: T) -> [T; 4] {
    let dx = x2 - x1;
    let divdif = (y2 - y1) / dx;
    let dzzdx = (divdif - s1) / dx;
    let dzdxx = (s2 - divdif) / dx;
    return [(dzdxx - dzzdx) / dx, dzzdx + dzzdx - dzdxx, s1, y1];
}

#[cfg(test)]
mod tests {
    use super::{fritschbutland_single_left, fritschbutland_single_middle, fritschbutland_single_right, slopes_fritschbutland};

    fn assert_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }

    #[test]
    fn single_intervals_match_full_slope_coefficients() {
        let xx = [0.0, 1.0, 2.0, 4.0];
        let yy = [0.0, 1.0, 1.5, 3.0];
        let slopes = slopes_fritschbutland(&xx, &yy);

        let left = fritschbutland_single_left(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2]);
        let middle = fritschbutland_single_middle(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);
        let right = fritschbutland_single_right(xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);

        assert_close(left[2], slopes[0]);
        assert_close(middle[2], slopes[1]);
        assert_close(right[2], slopes[2]);
    }
}
