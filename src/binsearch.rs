// binsearch.rs (splines library)
//! binary search algorithm to find the interval containing a search value; assumes monotonicity.

use num_traits::Float;

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/

/// For simplicity, does not accept NaNs or Inf values
/// Assumes the data are monotonous with the index.
/// | Argument | value |
/// |---|---|
/// | size | the number of intervals (nmax) |
/// | sval | search value |
/// | locator | a function returning  break values for a given index |
pub fn binary_search_interval<T: Float>(size: usize, sval: &T, locator: impl Fn(usize)->T)->Option<usize>{
	if _check_float(sval) {return None;}
	let mut n0 = 0;
	let mut n1 = size-1;
	let mut v0 = locator(n0);
	let v1 = locator(n1);
	if _check_float(&v0) || _check_float(&v1){return None;}
	if !interval_inside(sval, (&v0, &v1)){return None;}
	while (n1-n0) > 1 {
		let n = (n0+n1) / 2;
		let v = locator(n);
		if _check_float(&v){return None;}
		if interval_inside(sval, (&v0,&v)){
			n1 = n;
			//v1 = v;
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

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/

/// Returns `true` if a float value is neither NAN or INF.
fn _check_floats<T: Float, const N: usize>(vals: &[T;N])->bool{
	for val in vals.iter(){
		if val.is_nan() || val.is_infinite(){return true;}
	}
	return false;
}

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/

/// Iterate over val(k)-val(k-1) difference values in an iterator.
pub fn diff<T: Float>(slc: &[T])->impl Iterator<Item=T> + '_{
	return slc.windows(2).map(|w| w[1]-w[0]);
}

/// for a kernel (an array of fixed coefficients), a sum of products val(k)*coeff(k) + val(k-1)*coeff(k-1) + ... is calculated as an iterator over k. A lot of calculations using numerical methods can be formulated as a kernel transformation.
pub fn kernel_conv<'a, T: Float>(slc: &'a [T], kernel: &'a [T])->impl Iterator<Item=T> +'a {
	return slc.windows(kernel.len()).map(|w| _kernel_mult(w, kernel));
}

/// Apply kernel moving multiplication (folding).
fn _kernel_mult<T: Float>(window: &[T], kernel: &[T])->T{
	return window.iter().zip(kernel.iter()).map(|(w,k)| *w**k).fold(T::zero(), |acc, num| acc + num);
}
/// Check if a value is inside the interval.
pub fn interval_inside<T: Float>(val: &T, vals: (&T,&T))->bool{
	if val == vals.0 || val == vals.1 {return true;}
	let b1 = val > vals.0;
	let b2 = val > vals.1;
	return (vals.1 > vals.0 && b1 && !b2) || (vals.0 > vals.1 && b2 && !b1);
}

/// A simple implementation of a locator : locate a value inside a slice.
pub fn slice_locator<T: Float>(slc: &[T], loc: usize)->T {
	return slc[loc];
}
/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/