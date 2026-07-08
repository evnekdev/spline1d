// spline1d::spline.rs

//! The main structure in this module, PPData, encapsulates the calculated cubic spline coefficients and automatically handles interval location and interpolation using the appropriate set of cubic coefficients.

use std::fmt;
use std::collections::{HashMap};
use std::hash::{Hash};
use std::borrow::{Borrow};
use std::fs::{File};
use std::path::{Path};
use std::error::{Error};

use core::ops::{Add};
use num::{Float, Zero};
use csv::{Reader};
//use serde::{Serialize, Deserialize};
//use bincode::{Encode, Decode};

use crate::binsearch::{binary_search_interval,interval_inside,kernel_conv};
use crate::makima::{makima};
use crate::pchip::{pchip};
use crate::solve::{calculate_root};


/********************************************************************************************************************/
/// One-dimensional Piecewise-Polynomial
/// TODO make a zero-copy type where xx and yy are borrowed instead of cloning
#[derive(Debug, Clone)]
pub struct Spline<T: Float + fmt::Debug>{
	pub breaks_x: Vec<T>,
	pub breaks_y: Vec<T>,
	pub coeffs: Vec<[T;4]>,
	last: usize,
}
/*
impl<T: Float + fmt::Debug> num::Zero for Spline<T>{
	
	fn zero() -> Self {
		return Self {
			breaks_x: Vec::with_capacity(0),
			breaks_y: Vec::with_capacity(0),
			coeffs: Vec::with_capacity(0),
			last: 0,
		};
	}
	
	fn is_zero(&self) -> bool {
		return self.breaks_x.len() == 0 && self.coeffs.len() == 0;
	}
	
}

impl<T: Float + fmt::Debug> Add<Self> for Spline<T> {
	type Output = Spline<T>;
	fn add(self, rhs: Self) -> Spline<T> {
		return Self {
			breaks_x: Vec::new(),
			breaks_y: Vec::new(),
			coeffs: Vec::new(),
			last: 0,
		};
	}
	
}
*/
impl<T: Float + fmt::Debug> Spline<T> {
	
	pub fn new(xx: &[T], yy: &[T], ss: &[T])->Self { // xx is the principal variable, yy is a dependent variable
		let dxx : Vec<T> = kernel_conv(xx, &[-T::one(),T::one()]).collect(); // calculate differences in xx
		let dyy : Vec<T> = kernel_conv(yy, &[-T::one(),T::one()]).collect(); // calculate differences in yy
		let divdif : Vec<T> = dxx.iter().zip(dyy.iter()).map(|(x,y)| *y/ *x).collect(); // divide y by x
		//println!("divdif = {:?}", &divdif);
		let dzzdx : Vec<T> = dxx.iter().zip(divdif.iter().zip(ss.iter())).map(|(dx,(dydx,s))| (*dydx-*s)/ *dx).collect();
		let dzdxx : Vec<T> = dxx.iter().zip(divdif.iter().zip(ss[1..].iter())).map(|(dx,(dydx,s))| (*s-*dydx)/ *dx).collect();
		//println!("dzzdx = {:?}\ndzdxx = {:?}\n", &dzzdx, &dzdxx);
		let mut coeffs : Vec<[T;4]> = (0..xx.len()-1).map(|idx| [(dzdxx[idx]-dzzdx[idx])/dxx[idx],(dzzdx[idx]+dzzdx[idx])-dzdxx[idx],ss[idx],yy[idx]]).collect();
		//println!("coeffs = {:?}\n", &coeffs);
		for k in 0..coeffs.len(){
			for m in 0..4 {
				if coeffs[k][m].is_nan(){
					coeffs[k][m] = T::zero();
				}
			}
		}
		return Self{
			//breaks_x: xx.iter().map(|idx| *idx).collect(), // breaks at xx values
			//breaks_y: yy.iter().map(|idx| *idx).collect(),
			breaks_x: xx.to_vec(),
			breaks_y: yy.to_vec(),
			coeffs: coeffs,
			last: 0,
		};
	}
	/// Initializes a new instance using the makima method
	pub fn new_makima(xx: &[T], yy: &[T])->Self {
		return makima(xx, yy);
	}
	
	pub fn new_pchip(xx: &[T], yy: &[T])->Self {
		return pchip(xx, yy);
	}
	
	/// returns the index of the interval containing value x
	pub fn index(&self, x: &T)->Option<usize>{
		return binary_search_interval(self.breaks_x.len(), x, |loc| self.breaks_x[loc]);
	}
	/// checks whether value x is inside interval at idx
	fn check_index(&self, idx: &usize, x: &T)->bool{
		return interval_inside(x, (&self.breaks_x[*idx], &self.breaks_x[*idx+1]));
	}
	
	pub fn yvalue(&self, index: usize)-> T {
		return self.breaks_y[index];
	}
	
	pub fn max_value_index(&self)->usize {
		// linear scan of nodes
		let mut curr = self.coeffs[0][3];
		let mut currindex : usize = 0;
		for k in 0..self.breaks_y.len(){
			let val = self.breaks_y[k];
			if val > curr {
				curr = val;
				//println!("curr = {:?}", &curr);
				currindex = k;
			}
		}
		return currindex;
	}
	
	/// interpolate y for a given x
	pub fn interpolate(&self, x: &T)->Option<T>{
		let index = if self.check_index(&self.last, x){self.last} else {self.index(x)?};
		let xs = *x - self.breaks_x[index];
		return Some(self.coeffs[index].iter().fold(T::zero(), |acc, c| xs*acc + *c));
	}
	
	pub fn interpolate_linear(&self, x: &T)->Option<T>{
		let index = if self.check_index(&self.last, x){self.last} else {self.index(x)?};
		let x0 = self.breaks_x[index];
		let x1 = self.breaks_x[index+1];
		let xs = (*x - self.breaks_x[index])/(x1-x0);
		let y0 = self.breaks_y[index];
		let y1 = self.breaks_y[index+1];
		//println!("index = {:?}, x = {:?}, y0 = {:?}, y1 = {:?}", &index, &x, &y0, &y1);
		return Some(y0*(T::one()-xs) + y1*xs);
	}
	
	pub fn interpolate_diff1(&self, x: &T)->Option<T>{
		let index = if self.check_index(&self.last, x){self.last} else {self.index(x)?};
		return self.interpolate_diff1_for_index(x, index);
	}
	
	pub fn interpolate_diff2(&self, x: &T)->Option<T>{
		let index = if self.check_index(&self.last, x){self.last} else {self.index(x)?};
		return self.interpolate_diff2_for_index(x, index);
	}
	
	pub fn interpolate_for_index(&self, x: &T, index: usize)->Option<T>{
		let xs = *x - self.breaks_x[index];
		return Some(self.coeffs[index].iter().fold(T::zero(), |acc, c| xs*acc + *c));
	}
	
	pub fn interpolate_diff1_for_index(&self, x: &T, index: usize)->Option<T>{
		let xs = *x - self.breaks_x[index];
		let coeffs = &self.coeffs[index];
		let c1 = coeffs[0]*xs*xs;
		let c2 = coeffs[1]*xs;
		return Some(c1 + c1 + c1 + c2 + c2 + coeffs[2]);
	}
	
	pub fn interpolate_diff2_for_index(&self, x: &T, index: usize)->Option<T>{
		let xs = *x - self.breaks_x[index];
		let coeffs = &self.coeffs[index];
		let c1 = coeffs[0]*xs;
		let c2 = coeffs[1];
		return Some(c1 + c1 + c1 + c1 + c1 + c1 + c2 + c2);
	}
	
}

impl Spline<f64>{
	
	pub fn intersection_with(&self, other: &Spline<f64>, index1: usize, index2: usize)->Option<(f64,f64)>{
		let mut coeffs = vec![0.0;4];
		coeffs[0] = self.coeffs[index1][0]-other.coeffs[index2][0];
		coeffs[1] = self.coeffs[index1][1]-other.coeffs[index2][1];
		coeffs[2] = self.coeffs[index1][2]-other.coeffs[index2][2];
		coeffs[3] = self.coeffs[index1][3]-other.coeffs[index2][3];
		//println!("calculate_root = {:?}", &calculate_root(&coeffs));
		let x1 = self.breaks_x[index1] + calculate_root(&coeffs);
		let y1 = self.interpolate_for_index(&x1, index1)?;
		return Some((x1,y1));
	}
	
	pub fn intersection_with1(&self, other: &Spline<f64>, index1: usize, index2: usize)->Option<(f64,f64)>{
		//println!("intersection_with1");
		let mut x1 = self.breaks_x[index1];
		let mut x2 = other.breaks_x[index2];
		let func = |x| self.interpolate(&x).unwrap()-other.interpolate(&x).unwrap();
		let mut y1 = func(x1);
		let mut y2 = func(x2);
		while (x1-x2).abs() > 1.0e-6 {
			//println!("[{:?}-{:?}]", &x1, &x2);
			let x = (x1+x2)/2.0;
			let y = func(x);
			if ((y > 0.0) && (y1 > 0.0)) || ((y < 0.0) && (y1 < 0.0)){
				x1 = x;
			} else {
				x2 = x;
			}
		}
		let x = (x1+x2)/2.0;
		let y = self.interpolate(&x)?;
		return Some((x, y));
	}
	
}

