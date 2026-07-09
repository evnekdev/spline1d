// spline1d::methods::mod.rs
//! All methods of calculating cubic coefficients are grouped here.
//! 
//! Currently, the following local cubic interpolation methods are implemented:
//! `akima`, `makima`, `pchip`, `steffen`, `catmullrom`, `cardinal`, and
//! `fritschbutland`.

pub mod akima;
pub mod makima;
pub mod pchip;
pub mod steffen;
pub mod catmullrom;
pub mod cardinal;
pub mod fritschbutland;

pub use self::akima::{
    akima,
    slopes_akima,
    akima_single_left,
    akima_single_middle,
    akima_single_right,
    akima_single_left_alpha,
    akima_single_middle_alpha,
    akima_single_right_alpha,
};
pub use self::makima::{
    makima,
    slopes_makima,
    makima_single_left,
    makima_single_middle,
    makima_single_right,
    makima_single_left_alpha,
    makima_single_middle_alpha,
    makima_single_right_alpha,
};
pub use self::pchip::{
    pchip,
    slopes_pchip,
    pchip_single_left,
    pchip_single_middle,
    pchip_single_right,
    pchip_single_left_alpha,
    pchip_single_middle_alpha,
    pchip_single_right_alpha,
};
pub use self::steffen::{
    steffen,
    slopes_steffen,
    steffen_single_left,
    steffen_single_middle,
    steffen_single_right,
    steffen_single_left_alpha,
    steffen_single_middle_alpha,
    steffen_single_right_alpha,
};
pub use self::catmullrom::{
    catmullrom,
    slopes_catmullrom,
    catmullrom_single_left,
    catmullrom_single_middle,
    catmullrom_single_right,
    catmullrom_single_left_alpha,
    catmullrom_single_middle_alpha,
    catmullrom_single_right_alpha,
};
pub use self::cardinal::{
    cardinal,
    slopes_cardinal,
    cardinal_single_left,
    cardinal_single_middle,
    cardinal_single_right,
    cardinal_single_left_alpha,
    cardinal_single_middle_alpha,
    cardinal_single_right_alpha,
};
pub use self::fritschbutland::{
    fritschbutland,
    slopes_fritschbutland,
    fritschbutland_single_left,
    fritschbutland_single_middle,
    fritschbutland_single_right,
    fritschbutland_single_left_alpha,
    fritschbutland_single_middle_alpha,
    fritschbutland_single_right_alpha,
};
