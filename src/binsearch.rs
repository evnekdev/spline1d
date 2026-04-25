// binsearch.rs (splines library)

//! binary search algorithm to find the interval containing a search value; assumes monotonicity

use num::{Float};

use crate::{_interval_inside};

/************************************************************************************************************/

/// For simplicity, does not accept NaNs or Inf values
/// Assumes the data are monotonous with the index
pub fn binary_search_interval<T: Float>(size: usize, sval: &T, locator: impl Fn(usize)->T)->Option<usize>{
	if _check_float(sval) {return None;}
	let mut n0 = 0;
	let mut n1 = size-1;
	let mut v0 = locator(n0);
	let mut v1 = locator(n1);
	if _check_float(&v0) || _check_float(&v1){return None;}
	if !_interval_inside(sval, (&v0, &v1)){return None;}
	while (n1-n0) > 1 {
		let n = (n0+n1) / 2;
		let v = locator(n);
		if _check_float(&v){return None;}
		if _interval_inside(sval, (&v0,&v)){
			n1 = n;
			v1 = v;
		} else {
			n0 = n;
			v0 = v;
		}
	}
	return Some(n0);
}

fn _check_float<T: Float>(val: &T)->bool{
	return val.is_nan() || val.is_infinite();
}

/************************************************************************************************************/

/// TODO
pub fn binary_search_interval_nd<T: Float, const N: usize>(sizes: &[usize;N], svals: &[T;N], locator: impl Fn(&[usize;N])->[T;N])->Option<[usize;N]>{
	if _check_floats(svals){return None;}
	
	todo!();
}

fn _check_floats<T: Float, const N: usize>(vals: &[T;N])->bool{
	for val in vals.iter(){
		if val.is_nan() || val.is_infinite(){return true;}
	}
	return false;
}