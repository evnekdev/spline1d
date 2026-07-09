// catmullrom.rs (splines library)

//! Catmull-Rom cubic interpolation method.
//!
//! This is the zero-tension specialization of the cardinal method implemented
//! in this crate. Interior node derivatives are central differences:
//!
//! `m_i = (y_{i+1} - y_{i-1}) / (x_{i+1} - x_{i-1})`.

use num_traits::Float;

#[cfg(feature = "alloc")]
use alloc::vec;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

use crate::alpha::cubic_coeffs_to_alpha;
#[cfg(feature = "alloc")]
use crate::spline::Spline;

use super::cardinal::{
    cardinal_single_left,
    cardinal_single_middle,
    cardinal_single_right,
};
#[cfg(feature = "alloc")]
use super::cardinal::{cardinal, slopes_cardinal};

/// This function accepts x-values and y-values arrays and returns a spline interpolation container.
#[cfg(feature = "alloc")]
pub fn catmullrom<T: Float + core::fmt::Debug>(xx: &[T], yy: &[T]) -> Spline<T> {
    return cardinal(xx, yy, T::zero());
}

/// Estimation of tangent lines at `xx` points using Catmull-Rom central differences.
///
/// This is equivalent to [`slopes_cardinal`] with zero tension.
#[cfg(feature = "alloc")]
pub fn slopes_catmullrom<T: Float>(xx: &[T], yy: &[T]) -> Vec<T> {
    return slopes_cardinal(xx, yy, T::zero());
}

/// Cubic coefficients for a single left Catmull-Rom interval `[x1, x2]`.
///
/// Uses the first three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn catmullrom_single_left<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    return cardinal_single_left(x1, y1, x2, y2, x3, y3, T::zero());
}

/// Cubic coefficients for a single middle Catmull-Rom interval `[x1, x2]`.
///
/// Uses one neighbouring point on each side of the target interval. The returned
/// coefficients `[a, b, c, d]` are evaluated as `((a * dx + b) * dx + c) * dx + d`,
/// where `dx = x - x1`.
pub fn catmullrom_single_middle<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    return cardinal_single_middle(x0, y0, x1, y1, x2, y2, x3, y3, T::zero());
}

/// Cubic coefficients for a single right Catmull-Rom interval `[x1, x2]`.
///
/// Uses the last three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn catmullrom_single_right<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 4] {
    return cardinal_single_right(x0, y0, x1, y1, x2, y2, T::zero());
}

/// Normalized `[alpha0, alpha1]` coefficients for a single left Catmull-Rom interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn catmullrom_single_left_alpha<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = catmullrom_single_left(x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single middle Catmull-Rom interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn catmullrom_single_middle_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = catmullrom_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single right Catmull-Rom interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn catmullrom_single_right_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 2] {
    let coeffs = catmullrom_single_right(x0, y0, x1, y1, x2, y2);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

#[cfg(test)]
mod tests {
    use super::{catmullrom_single_left, catmullrom_single_middle, catmullrom_single_right, slopes_catmullrom};

    fn assert_close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }

    #[test]
    fn single_intervals_match_full_slope_coefficients() {
        let xx = [0.0, 1.0, 2.0, 4.0];
        let yy = [0.0, 1.0, 1.5, 3.0];
        let slopes = slopes_catmullrom(&xx, &yy);

        let left = catmullrom_single_left(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2]);
        let middle = catmullrom_single_middle(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);
        let right = catmullrom_single_right(xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);

        assert_close(left[2], slopes[0]);
        assert_close(middle[2], slopes[1]);
        assert_close(right[2], slopes[2]);
    }
}
