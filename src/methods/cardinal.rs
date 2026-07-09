// cardinal.rs (splines library)

//! Cardinal spline interpolation method.
//!
//! This is the tension-controlled generalization of Catmull-Rom.  The parameter
//! `tension` is used in the common form
//!
//! `d[i] = (1 - tension) * (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])`.
//!
//! Therefore `tension = 0` gives Catmull-Rom, and increasing `tension` shrinks
//! the node derivatives.

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
///
/// `tension = 0` gives Catmull-Rom.  Larger tension reduces all tangents by the
/// factor `1 - tension`.
#[cfg(feature = "alloc")]
pub fn cardinal<T: Float + core::fmt::Debug>(xx: &[T], yy: &[T], tension: T) -> Spline<T> {
    let ss = slopes_cardinal(xx, yy, tension);
    return Spline::new(xx, yy, &ss);
}

/// Estimation of the tangent lines at `xx` points using Cardinal spline slopes.
///
/// Endpoint slopes are the adjacent secant slopes scaled by `1 - tension`.
/// Interior slopes use the two-sided neighbour slope scaled by `1 - tension`.
#[cfg(feature = "alloc")]
pub fn slopes_cardinal<T: Float>(xx: &[T], yy: &[T], tension: T) -> Vec<T> {
    let n = xx.len();

    assert!(n >= 2, "cardinal requires at least two points");
    assert!(yy.len() == n, "xx and yy must have equal length");

    let factor = T::one() - tension;
    let h: Vec<T> = diff(xx).collect();
    let delta: Vec<T> = diff(yy).zip(h.iter()).map(|(dy, dx)| dy / *dx).collect();

    if n == 2 {
        return vec![factor * delta[0]; n];
    }

    let mut slopes = vec![T::zero(); n];
    slopes[0] = factor * delta[0];
    slopes[n - 1] = factor * delta[n - 2];

    for i in 1..(n - 1) {
        slopes[i] = factor * (yy[i + 1] - yy[i - 1]) / (xx[i + 1] - xx[i - 1]);
    }

    return slopes;
}

/// Cubic coefficients for a single left Cardinal interval `[x1, x2]`.
///
/// Uses the first three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn cardinal_single_left<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, tension: T) -> [T; 4] {
    let factor = T::one() - tension;
    let s1 = factor * secant(x1, y1, x2, y2);
    let s2 = factor * secant(x1, y1, x3, y3);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single middle Cardinal interval `[x1, x2]`.
///
/// Uses one neighbouring point on each side of the target interval. The returned
/// coefficients `[a, b, c, d]` are evaluated as `((a * dx + b) * dx + c) * dx + d`,
/// where `dx = x - x1`.
pub fn cardinal_single_middle<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, tension: T) -> [T; 4] {
    let factor = T::one() - tension;
    let s1 = factor * secant(x0, y0, x2, y2);
    let s2 = factor * secant(x1, y1, x3, y3);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single right Cardinal interval `[x1, x2]`.
///
/// Uses the last three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn cardinal_single_right<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, tension: T) -> [T; 4] {
    let factor = T::one() - tension;
    let s1 = factor * secant(x0, y0, x2, y2);
    let s2 = factor * secant(x1, y1, x2, y2);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single left Cardinal interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn cardinal_single_left_alpha<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, tension: T) -> [T; 2] {
    let coeffs = cardinal_single_left(x1, y1, x2, y2, x3, y3, tension);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single middle Cardinal interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn cardinal_single_middle_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, tension: T) -> [T; 2] {
    let coeffs = cardinal_single_middle(x0, y0, x1, y1, x2, y2, x3, y3, tension);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single right Cardinal interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn cardinal_single_right_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, tension: T) -> [T; 2] {
    let coeffs = cardinal_single_right(x0, y0, x1, y1, x2, y2, tension);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

#[inline]
fn secant<T: Float>(xa: T, ya: T, xb: T, yb: T) -> T {
    return (yb - ya) / (xb - xa);
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
    use super::{cardinal_single_left, cardinal_single_middle, cardinal_single_right, slopes_cardinal};

    fn assert_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }

    #[test]
    fn single_intervals_match_full_slope_coefficients() {
        let xx = [0.0, 1.0, 2.0, 4.0];
        let yy = [0.0, 1.0, 1.5, 3.0];
        let tension = 0.25;
        let slopes = slopes_cardinal(&xx, &yy, tension);

        let left = cardinal_single_left(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], tension);
        let middle = cardinal_single_middle(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], xx[3], yy[3], tension);
        let right = cardinal_single_right(xx[1], yy[1], xx[2], yy[2], xx[3], yy[3], tension);

        assert_close(left[2], slopes[0]);
        assert_close(middle[2], slopes[1]);
        assert_close(right[2], slopes[2]);
    }
}
