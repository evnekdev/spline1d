// makima.rs (splines library)

//! Makima interpolation method (copied from Matlab)

use num_traits::Float;

#[cfg(feature = "alloc")]
use alloc::vec;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

#[cfg(feature = "alloc")]
use crate::spline::Spline;
use crate::binsearch::{diff, kernel_conv};
use crate::alpha::cubic_coeffs_to_alpha;

/// This function accepts x-values and y-values arrays and returns a spline interpolation container
#[cfg(feature = "alloc")]
pub fn makima<T: Float + core::fmt::Debug>(xx: &[T], yy: &[T])->Spline<T>{
	let ss = slopes_makima(xx, yy);
	return Spline::new(xx, yy, &ss);
}

/// Estimation of the tangent lines at xx points using the makima method
#[cfg(feature = "alloc")]
pub fn slopes_makima<T: Float + core::fmt::Debug>(xx: &[T], yy: &[T])->Vec<T>{
	let hh : Vec<T> = diff(xx).collect();
	let delta : Vec<T> = diff(yy).zip(hh.iter()).map(|(dy,dx)| dy/ *dx).collect();
	// special case of two points, use linear slope
	if xx.len() == 2 {
		return vec![delta[0];2];
	}
	// calculate the missing deltas for points 0, 1, n, n-1.
	let n = xx.len();
	let delta_m1 : T = (delta[0]+delta[0]) - delta[1];
	let delta_m2 : T = (delta_m1+delta_m1) - delta[0];
	let delta_n  : T = (delta[n-2]+delta[n-2]) - delta[n-3];
	let delta_n1 : T = (delta_n+delta_n) - delta[n-2];
	let delta_prefix = vec![delta_m2,delta_m1];
	let delta_suffix = vec![delta_n, delta_n1];
	let delta_new : Vec<T> = delta_prefix.iter().chain(delta.iter()).chain(delta_suffix.iter()).map(|&v| v).collect();
	let k1 = [-T::one(),T::one()];
	let half : T = T::one()/(T::one() + T::one());
	let k2 = [half, half];
	let it1 = kernel_conv(&delta_new, &k1).map(|v| v.abs());
	let it2 = kernel_conv(&delta_new, &k2).map(|v| v.abs());
	let weights : Vec<T> = it1.zip(it2).map(|(v1,v2)| v1+v2).collect();
	let k3 = [T::one(), T::zero(), T::one()];
	let weights12 : Vec<T> = kernel_conv(&weights, &k3).collect();
	let s1 : Vec<T> = weights[2..].iter().zip(delta_new[1..n+1].iter()).map(|(w,d)| *w* *d).collect();
	let s2 : Vec<T> = weights[0..n].iter().zip(delta_new[2..n+2].iter()).map(|(w,d)| *w* *d).collect();
	let ss : Vec<T> = weights12.iter().zip(s1.iter().zip(s2.iter())).map(|(w,(s1,s2))| (*s1+*s2)/ *w).collect();
	return ss;
}

/// Cubic coefficients for a single left Makima interval `[x1, x2]`.
///
/// Uses only the first three data points and Makima's linear extrapolation of
/// missing secant slopes on the left. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn makima_single_left<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let d12 = secant(x1, y1, x2, y2);
    let d23 = secant(x2, y2, x3, y3);

    let d0 = (d12 + d12) - d23;
    let dm1 = (d0 + d0) - d12;
    let d34 = (d23 + d23) - d12;

    let s1 = makima_slope(dm1, d0, d12, d23);
    let s2 = makima_slope(d0, d12, d23, d34);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single middle Makima interval `[x1, x2]`.
///
/// Uses one neighbouring point on each side of the target interval. The returned
/// coefficients `[a, b, c, d]` are evaluated as `((a * dx + b) * dx + c) * dx + d`,
/// where `dx = x - x1`.
pub fn makima_single_middle<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 4] {
    let d01 = secant(x0, y0, x1, y1);
    let d12 = secant(x1, y1, x2, y2);
    let d23 = secant(x2, y2, x3, y3);

    let dm1 = (d01 + d01) - d12;
    let d34 = (d23 + d23) - d12;

    let s1 = makima_slope(dm1, d01, d12, d23);
    let s2 = makima_slope(d01, d12, d23, d34);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Cubic coefficients for a single right Makima interval `[x1, x2]`.
///
/// Uses only the last three data points and Makima's linear extrapolation of
/// missing secant slopes on the right. The returned coefficients `[a, b, c, d]`
/// are evaluated as `((a * dx + b) * dx + c) * dx + d`, where `dx = x - x1`.
pub fn makima_single_right<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 4] {
    let d01 = secant(x0, y0, x1, y1);
    let d12 = secant(x1, y1, x2, y2);

    let dm1 = (d01 + d01) - d12;
    let d23 = (d12 + d12) - d01;
    let d34 = (d23 + d23) - d12;

    let s1 = makima_slope(dm1, d01, d12, d23);
    let s2 = makima_slope(d01, d12, d23, d34);

    return cubic_coeffs(x1, y1, x2, y2, s1, s2);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single left Makima interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn makima_single_left_alpha<T: Float>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = makima_single_left(x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single middle Makima interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn makima_single_middle_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T) -> [T; 2] {
    let coeffs = makima_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

/// Normalized `[alpha0, alpha1]` coefficients for a single right Makima interval `[x1, x2]`.
///
/// The returned coefficients are for
/// `y = y1*(1-t) + y2*t + (1-t)*t*(alpha0 + alpha1*t)`,
/// where `t = (x - x1) / (x2 - x1)`.
pub fn makima_single_right_alpha<T: Float>(x0: T, y0: T, x1: T, y1: T, x2: T, y2: T) -> [T; 2] {
    let coeffs = makima_single_right(x0, y0, x1, y1, x2, y2);
    return cubic_coeffs_to_alpha(coeffs, x2 - x1);
}

#[inline]
fn secant<T: Float>(xa: T, ya: T, xb: T, yb: T) -> T {
    return (yb - ya) / (xb - xa);
}

#[inline]
fn makima_slope<T: Float>(dm2: T, dm1: T, d0: T, d1: T) -> T {
    let half = T::one() / (T::one() + T::one());
    let w1 = (d1 - d0).abs() + (half * (d1 + d0)).abs();
    let w2 = (dm1 - dm2).abs() + (half * (dm1 + dm2)).abs();
    let w = w1 + w2;

    if w == T::zero() {
        return (d0 + dm1) * half;
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
