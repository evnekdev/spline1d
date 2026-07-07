// lib.rs (splines library)

//! This crate performs cubic spline interpolation in pure Rust (1-D curves); currently, two methods of cubic interpolation are available: makima (fully developed and tested) and pchip (in progress). The original publications of the cubic interpolation methods are located in the project folder; the working formulas for makima are copied after the Matlab in-built methods.
//!
//! The root contains some common utility functions used elsewhere


pub mod pchip;
pub mod makima;
pub mod binsearch;
pub mod spline;
pub mod multispline;
pub mod solve;
pub mod search_tree;

use num::{Float};

pub use crate::makima::{makima};
pub use crate::pchip::{pchip};
pub use crate::binsearch::{binary_search_interval};
pub use crate::spline::{Spline};
pub use crate::multispline::{MultiSpline, load_mpp_from_csv};
pub use crate::search_tree::{SearchNode,SearchTree};

/// Iterate over val(k)-val(k-1) difference values in an iterator
pub fn diff<T: Float>(slc: &[T])->impl Iterator<Item=T> + '_{
	return slc.windows(2).map(|w| w[1]-w[0]);
}

/// for a kernel (an array of fixed coefficients), a sum of products val(k)*coeff(k) + val(k-1)*coeff(k-1) + ... is calculated as an iterator over k. A lot of calculations using numerical methods can be formulated as a kernel transformation
pub fn kernel_conv<'a, T: Float>(slc: &'a [T], kernel: &'a [T])->impl Iterator<Item=T> +'a {
	return slc.windows(kernel.len()).map(|w| _kernel_mult(w, kernel));
}


fn _kernel_mult<T: Float>(window: &[T], kernel: &[T])->T{
	return window.iter().zip(kernel.iter()).map(|(w,k)| *w**k).fold(T::zero(), |acc, num| acc + num);
}
/// check if a value is inside the interval
fn _interval_inside<T: Float>(val: &T, vals: (&T,&T))->bool{
	if val == vals.0 || val == vals.1 {return true;}
	let b1 = val > vals.0;
	let b2 = val > vals.1;
	return (vals.1 > vals.0 && b1 && !b2) || (vals.0 > vals.1 && b2 && !b1);
}

/// Locate a value inside a slice
pub fn slice_locator<T: Float>(slc: &[T], loc: usize)->T {
	return slc[loc];
}

/********************************************************************************************************************/
