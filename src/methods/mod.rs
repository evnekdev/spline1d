// spline1d::methods::mod.rs

pub mod makima;
pub mod pchip;
pub mod steffen;


pub use crate::makima::{makima, makima_single_left, makima_single_middle, makima_single_right, makima_single_left_alpha, makima_single_middle_alpha, makima_single_right_alpha};
pub use crate::pchip::{pchip, pchip_single_left, pchip_single_middle, pchip_single_right, pchip_single_left_alpha, pchip_single_middle_alpha, pchip_single_right_alpha};
pub use crate::steffen::{steffen, slopes_steffen, steffen_single_left, steffen_single_middle, steffen_single_right, steffen_single_left_alpha, steffen_single_middle_alpha, steffen_single_right_alpha};