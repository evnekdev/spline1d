// pchip.rs (splines library)

//! PCHIP interpolation method (translated from Matlab's `pchipSlopes`).

use num::Float;

use crate::binsearch::diff;
use crate::spline::Spline;
use crate::alpha::cubic_coeffs_to_alpha;

/// This function accepts x-values and y-values arrays and returns a spline interpolation container.
pub fn pchip<T: Float + std::fmt::Debug>(xx: &[T], yy: &[T]) -> Spline<T> {
    let hh: Vec<T> = diff(xx).collect();
    let delta: Vec<T> = diff(yy).zip(hh.iter()).map(|(dy, dx)| dy / *dx).collect();
    let ss = slopes_pchip(xx, yy, &delta);
    return Spline::new(xx, yy, &ss);
}

/// Estimation of the tangent lines at `xx` points using the PCHIP method.
///
/// This is a direct real-valued translation of Matlab's internal `pchipSlopes`
/// routine. For a monotone sequence, the slopes are shape-preserving: interior
/// slopes are zero when adjacent secant slopes change sign, and otherwise use
/// the weighted harmonic mean formula.
pub fn slopes_pchip<T: Float>(xx: &[T], _yy: &[T], delta: &[T]) -> Vec<T> {
    let n = xx.len();

    assert!(n >= 2, "pchip requires at least two points");
    assert!(delta.len() == n - 1, "delta must have length xx.len() - 1");

    // Special case n = 2, use linear interpolation.
    if n == 2 {
        return vec![delta[0]; n];
    }

    let zero = T::zero();
    let one = T::one();
    let two = one + one;
    let three = two + one;

    let h: Vec<T> = diff(xx).collect();
    let mut d = vec![zero; n];

    // Slopes at interior points.
    // d(k) = weighted average of delta(k-1) and delta(k) when they have the
    // same sign; d(k) = 0 when they have opposite signs or either is zero.
    for i in 0..(n - 2) {
        let del1 = delta[i];
        let del2 = delta[i + 1];

        if del1 * del2 > zero {
            let hs = h[i] + h[i + 1];
            let w1 = (h[i] + hs) / (three * hs);
            let w2 = (hs + h[i + 1]) / (three * hs);

            let abs1 = del1.abs();
            let abs2 = del2.abs();
            let dmax = if abs1 > abs2 { abs1 } else { abs2 };
            let dmin = if abs1 < abs2 { abs1 } else { abs2 };

            d[i + 1] = dmin / (w1 * (del1 / dmax) + w2 * (del2 / dmax));
        }
    }

    // Slopes at end points.
    // Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.
    d[0] = ((two * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);
    if d[0] * delta[0] <= zero {
        d[0] = zero;
    } else if delta[0] * delta[1] < zero && d[0].abs() > (three * delta[0]).abs() {
        d[0] = three * delta[0];
    }

    d[n - 1] = ((two * h[n - 2] + h[n - 3]) * delta[n - 2] - h[n - 2] * delta[n - 3])
        / (h[n - 2] + h[n - 3]);
    if d[n - 1] * delta[n - 2] <= zero {
        d[n - 1] = zero;
    } else if delta[n - 2] * delta[n - 3] < zero && d[n - 1].abs() > (three * delta[n - 2]).abs() {
        d[n - 1] = three * delta[n - 2];
    }

    return d;
}

/// Cubic coefficients for a single left PCHIP interval `[x1, x2]`.
///
/// Uses the first three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn pchip_single_left<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let h12 = x2 - x1;
    let h23 = x3 - x2;
    let d12 = (y2 - y1) / h12;
    let d23 = (y3 - y2) / h23;

    let s1 = pchip_endpoint_slope(h12, h23, d12, d23);
    let s2 = pchip_interior_slope(h12, h23, d12, d23);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single middle PCHIP interval `[x1, x2]`.
///
/// Uses one neighbouring point on each side of the target interval. The returned
/// coefficients `[a, b, c, d]` are evaluated as `((a * dx + b) * dx + c) * dx + d`,
/// where `dx = x - x1`.
pub fn pchip_single_middle<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let h01 = x1 - x0;
    let h12 = x2 - x1;
    let h23 = x3 - x2;
    let d01 = (y1 - y0) / h01;
    let d12 = (y2 - y1) / h12;
    let d23 = (y3 - y2) / h23;

    let s1 = pchip_interior_slope(h01, h12, d01, d12);
    let s2 = pchip_interior_slope(h12, h23, d12, d23);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single right PCHIP interval `[x1, x2]`.
///
/// Uses the last three data points. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn pchip_single_right<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 4] {
    let h01 = x1 - x0;
    let h12 = x2 - x1;
    let d01 = (y1 - y0) / h01;
    let d12 = (y2 - y1) / h12;

    let s1 = pchip_interior_slope(h01, h12, d01, d12);
    let s2 = pchip_endpoint_slope(h12, h01, d12, d01);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single left PCHIP interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn pchip_single_left_alpha<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = pchip_single_left(x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single middle PCHIP interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn pchip_single_middle_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = pchip_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single right PCHIP interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn pchip_single_right_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 2] {
    let coeffs = pchip_single_right(x0, y0, x1, y1, x2, y2);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

#[inline]
fn pchip_interior_slope<T: Float>(h1: T, h2: T, del1: T, del2: T) -> T {
    let zero = T::zero();

    if del1 * del2 <= zero {
        return zero;
    }

    let one = T::one();
    let three = one + one + one;
    let hs = h1 + h2;
    let w1 = (h1 + hs) / (three * hs);
    let w2 = (hs + h2) / (three * hs);
    let abs1 = del1.abs();
    let abs2 = del2.abs();
    let dmax = if abs1 > abs2 { abs1 } else { abs2 };
    let dmin = if abs1 < abs2 { abs1 } else { abs2 };

    return dmin / (w1 * (del1 / dmax) + w2 * (del2 / dmax));
}

#[inline]
fn pchip_endpoint_slope<T: Float>(h1: T, h2: T, del1: T, del2: T) -> T {
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
fn cubic_coeffs<T: Float>(x1: T, y1: T, x2: T, y2: T, s1: T, s2: T) -> [T; 4] {
    let dx = x2 - x1;
    let divdif = (y2 - y1) / dx;
    let dzzdx = (divdif - s1) / dx;
    let dzdxx = (s2 - divdif) / dx;
    return [(dzdxx - dzzdx) / dx, dzzdx + dzzdx - dzdxx, s1, y1];
}
