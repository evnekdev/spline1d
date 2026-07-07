// makima.rs
use std::mem::{size_of_val};
use spline1d::*;

pub fn main(){
	let xx = vec![0.0, 0.2, 0.4, 0.55, 0.65];
	let yy = vec![2845.2, 2688.0, 2448.2, 2114.4, 1807.8];
	let pp = makima(&xx, &yy);
	let xx : Vec<f64> = vec![0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55];
	for x in xx.iter(){
		let y = pp.interpolate(x);
		println!("{:?},{:?}", x, &y);
	}
}