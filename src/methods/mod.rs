// spline1d::methods::mod.rs
//! All methods of calculating cubic coefficients are grouped here.
//! 
//! Currently, the following local cubic interpolation methods are implemented:
//! `akima`, `makima`, `pchip`, `steffen`, `catmullrom`, `cardinal`, and
//! `fritschbutland`.

use num::{Float};

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

#[derive(Clone,Copy,Debug)]
pub enum InterpolationType<T: Float> {
	AKIMA,
	MAKIMA,
	PCHIP,
	STEFFEN,
	CATMULLROM,
	CARDINAL(T),
	FRITSCHBUTLAND,
}

/*******************************************************************************************************************************************************************/
/*******************************************************************************************************************************************************************/

/// Generic method to calculate cubic coefficients on a single left interval.
pub fn cubic_single_left<T: Float>(itype: InterpolationType<T>, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T)->[T;4]{
	match itype {
		InterpolationType::AKIMA             => {return akima_single_left(x1, y1, x2, y2, x3, y3);}
		InterpolationType::MAKIMA            => {return makima_single_left(x1, y1, x2, y2, x3, y3);}
		InterpolationType::PCHIP             => {return pchip_single_left(x1, y1, x2, y2, x3, y3);}
		InterpolationType::STEFFEN           => {return steffen_single_left(x1, y1, x2, y2, x3, y3);}
		InterpolationType::CATMULLROM        => {return catmullrom_single_left(x1, y1, x2, y2, x3, y3);}
		InterpolationType::CARDINAL(tension) => {return cardinal_single_left(x1, y1, x2, y2, x3, y3, tension);}
		InterpolationType::FRITSCHBUTLAND    => {return fritschbutland_single_left(x1, y1, x2, y2, x3, y3);}
	}
}

/// Generic method to calculate cubic coefficients on a single middle interval.
pub fn cubic_single_middle<T: Float>(itype: InterpolationType<T>, x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T)->[T;4] {
	match itype {
		InterpolationType::AKIMA             => {return akima_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::MAKIMA            => {return makima_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::PCHIP             => {return pchip_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::STEFFEN           => {return steffen_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::CATMULLROM        => {return catmullrom_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::CARDINAL(tension) => {return cardinal_single_middle(x0, y0, x1, y1, x2, y2, x3, y3, tension);}
		InterpolationType::FRITSCHBUTLAND    => {return fritschbutland_single_middle(x0, y0, x1, y1, x2, y2, x3, y3);}
	}
}

/// Generic method to calculate cubic coefficients on a single right interval.
pub fn cubic_single_right<T: Float>(itype: InterpolationType<T>, x0: T, y0: T, x1: T, y1: T, x2: T, y2: T)->[T;4] {
	match itype {
		InterpolationType::AKIMA             => {return akima_single_right(x0, y0, x1, y1, x2, y2);}
		InterpolationType::MAKIMA            => {return makima_single_right(x0, y0, x1, y1, x2, y2);}
		InterpolationType::PCHIP             => {return pchip_single_right(x0, y0, x1, y1, x2, y2);}
		InterpolationType::STEFFEN           => {return steffen_single_right(x0, y0, x1, y1, x2, y2);}
		InterpolationType::CATMULLROM        => {return catmullrom_single_right(x0, y0, x1, y1, x2, y2);}
		InterpolationType::CARDINAL(tension) => {return cardinal_single_right(x0, y0, x1, y1, x2, y2, tension);}
		InterpolationType::FRITSCHBUTLAND    => {return fritschbutland_single_right(x0, y0, x1, y1, x2, y2);}
	}
}

/// TODO
pub fn cubic_single_left_alpha<T: Float>(itype: InterpolationType<T>, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T)->[T;2]{
	match itype {
		InterpolationType::AKIMA             => {return akima_single_left_alpha(x1, y1, x2, y2, x3, y3);}
		InterpolationType::MAKIMA            => {return makima_single_left_alpha(x1, y1, x2, y2, x3, y3);}
		InterpolationType::PCHIP             => {return pchip_single_left_alpha(x1, y1, x2, y2, x3, y3);}
		InterpolationType::STEFFEN           => {return steffen_single_left_alpha(x1, y1, x2, y2, x3, y3);}
		InterpolationType::CATMULLROM        => {return catmullrom_single_left_alpha(x1, y1, x2, y2, x3, y3);}
		InterpolationType::CARDINAL(tension) => {return cardinal_single_left_alpha(x1, y1, x2, y2, x3, y3, tension);}
		InterpolationType::FRITSCHBUTLAND    => {return fritschbutland_single_left_alpha(x1, y1, x2, y2, x3, y3);}
	}
}

/// TODO
pub fn cubic_single_middle_alpha<T: Float>(itype: InterpolationType<T>, x0: T, y0: T, x1: T, y1: T, x2: T, y2: T, x3: T, y3: T)->[T;2] {
	match itype {
		InterpolationType::AKIMA             => {return akima_single_middle_alpha(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::MAKIMA            => {return makima_single_middle_alpha(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::PCHIP             => {return pchip_single_middle_alpha(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::STEFFEN           => {return steffen_single_middle_alpha(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::CATMULLROM        => {return catmullrom_single_middle_alpha(x0, y0, x1, y1, x2, y2, x3, y3);}
		InterpolationType::CARDINAL(tension) => {return cardinal_single_middle_alpha(x0, y0, x1, y1, x2, y2, x3, y3, tension);}
		InterpolationType::FRITSCHBUTLAND    => {return fritschbutland_single_middle_alpha(x0, y0, x1, y1, x2, y2, x3, y3);}
	}
}

/// TODO
pub fn cubic_single_right_alpha<T: Float>(itype: InterpolationType<T>, x0: T, y0: T, x1: T, y1: T, x2: T, y2: T)->[T;2] {
	match itype {
		InterpolationType::AKIMA             => {return akima_single_right_alpha(x0, y0, x1, y1, x2, y2);}
		InterpolationType::MAKIMA            => {return makima_single_right_alpha(x0, y0, x1, y1, x2, y2);}
		InterpolationType::PCHIP             => {return pchip_single_right_alpha(x0, y0, x1, y1, x2, y2);}
		InterpolationType::STEFFEN           => {return steffen_single_right_alpha(x0, y0, x1, y1, x2, y2);}
		InterpolationType::CATMULLROM        => {return catmullrom_single_right_alpha(x0, y0, x1, y1, x2, y2);}
		InterpolationType::CARDINAL(tension) => {return cardinal_single_right_alpha(x0, y0, x1, y1, x2, y2, tension);}
		InterpolationType::FRITSCHBUTLAND    => {return fritschbutland_single_right_alpha(x0, y0, x1, y1, x2, y2);}
	}
}

/*******************************************************************************************************************************************************************/
/*******************************************************************************************************************************************************************/