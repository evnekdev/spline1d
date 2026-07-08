// lib.rs (splines library)

//! This crate performs cubic spline interpolation in pure Rust (1-D curves); currently, three methods of cubic interpolation are available: makima, pchip, and steffen. The original publications of the cubic interpolation methods are located in the project folder; the working formulas for makima are copied after the Matlab in-built methods.
//!
//! The root contains some common utility functions used elsewhere


pub mod alpha;
pub mod pchip;
pub mod makima;
pub mod steffen;
pub mod binsearch;
pub mod spline;
pub mod multispline;
pub mod solve;
pub mod searchtree;

use num::{Float};

pub use crate::makima::{makima, makima_single_left, makima_single_middle, makima_single_right, makima_single_left_alpha, makima_single_middle_alpha, makima_single_right_alpha};
pub use crate::pchip::{pchip, pchip_single_left, pchip_single_middle, pchip_single_right, pchip_single_left_alpha, pchip_single_middle_alpha, pchip_single_right_alpha};
pub use crate::steffen::{steffen, slopes_steffen, steffen_single_left, steffen_single_middle, steffen_single_right, steffen_single_left_alpha, steffen_single_middle_alpha, steffen_single_right_alpha};
pub use crate::alpha::{cubic_coeffs_to_alpha, cubic_coeffs_to_alpha_unit, alpha_to_cubic_coeffs, alpha_to_standard_cubic_coeffs};
pub use crate::binsearch::{binary_search_interval};
pub use crate::spline::{Spline};
pub use crate::multispline::{MultiSpline, load_multispline_from_csv};
pub use crate::searchtree::{SearchNode,SearchTree};

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/


