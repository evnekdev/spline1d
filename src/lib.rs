// lib.rs (splines library)

//! This crate performs cubic spline interpolation in pure Rust (1-D curves).
//! Currently, three methods of cubic interpolation are available: `makima`, `pchip`, and `steffen`.
//! The original publications of the cubic interpolation methods are located in the project folder; the working formulas for makima and pchip are copied after the Matlab in-built methods.
//! 
//! Cubic interpolation is available in two flavors: single-interval functions calculating cubic coefficients directly (4-values) or producing a structure with interval breaks alongside with cubic coefficients used to lookup the containing interpolation intervals for `xs`.
//! Single-interval functions can produce either cubic coefficients (4-values) or normalized endpoints coefficients (alpha coefficients, 2-values).
//! 
//! An advanced multivariable 1D interpolation is available; interpolation structures for n variables x1, x2, ..., xn running along x0 can be calculated simultaneously and used to do any pair-wise lookup xi->xj.
//!
//! In case of non-monotonous variables, a reverse interpolation search exists: a `SearchTree` is produced to lookup all x values for a given y.

pub mod alpha;
pub mod methods;
pub mod binsearch;
pub mod spline;
pub mod multispline;
pub mod solve;
pub mod searchtree;

use num::{Float};

pub use crate::methods::{makima, makima_single_left, makima_single_middle, makima_single_right, makima_single_left_alpha, makima_single_middle_alpha, makima_single_right_alpha};
pub use crate::methods::{pchip, pchip_single_left, pchip_single_middle, pchip_single_right, pchip_single_left_alpha, pchip_single_middle_alpha, pchip_single_right_alpha};
pub use crate::methods::{steffen, slopes_steffen, steffen_single_left, steffen_single_middle, steffen_single_right, steffen_single_left_alpha, steffen_single_middle_alpha, steffen_single_right_alpha};
pub use crate::alpha::{cubic_coeffs_to_alpha, cubic_coeffs_to_alpha_unit, alpha_to_cubic_coeffs, alpha_to_standard_cubic_coeffs};
pub use crate::binsearch::{binary_search_interval};
pub use crate::spline::{Spline};
pub use crate::multispline::{MultiSpline, load_multispline_from_csv};
pub use crate::searchtree::{SearchNode,SearchTree};

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/


