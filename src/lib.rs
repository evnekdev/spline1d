// lib.rs (splines library)

//! This crate performs cubic spline interpolation in pure Rust (1-D curves); currently, two methods of cubic interpolation are available: makima and pchip. The original publications of the cubic interpolation methods are located in the project folder; the working formulas for makima are copied after the Matlab in-built methods.
//!
//! The root contains some common utility functions used elsewhere


pub mod pchip;
pub mod makima;
pub mod binsearch;
pub mod spline;
pub mod multispline;
pub mod solve;
pub mod searchtree;

use num::{Float};

pub use crate::makima::{makima, makima_single_left, makima_single_middle, makima_single_right};
pub use crate::pchip::{pchip, pchip_single_left, pchip_single_middle, pchip_single_right};
pub use crate::binsearch::{binary_search_interval};
pub use crate::spline::{Spline};
pub use crate::multispline::{MultiSpline, load_multispline_from_csv};
pub use crate::searchtree::{SearchNode,SearchTree};

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/


