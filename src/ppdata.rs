// ppdata.rs (splines library)

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
use serde::{Serialize, Deserialize};
use bincode::{Encode, Decode};

use crate::{kernel_conv};
use crate::{_interval_inside};
use crate::binsearch::{binary_search_interval};
use crate::makima::{makima};
use crate::solve::{calculate_root};

/********************************************************************************************************************/
/// One-dimensional Piecewise-Polynomial
/// TODO make a zero-copy type where xx and yy are borrowed instead of cloning
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct PPData<T: Float + fmt::Debug>{
	pub breaks_x: Vec<T>,
	pub breaks_y: Vec<T>,
	pub coeffs: Vec<[T;4]>,
	last: usize,
}

impl<T: Float + fmt::Debug> num::Zero for PPData<T>{
	
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

impl<T: Float + fmt::Debug> Add<Self> for PPData<T> {
	type Output = PPData<T>;
	fn add(self, rhs: Self) -> PPData<T> {
		return Self {
			breaks_x: Vec::new(),
			breaks_y: Vec::new(),
			coeffs: Vec::new(),
			last: 0,
		};
	}
	
}

impl<T: Float + fmt::Debug> PPData<T> {
	
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
			breaks_x: xx.iter().map(|idx| *idx).collect(), // breaks at xx values
			breaks_y: yy.iter().map(|idx| *idx).collect(),
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
		return _interval_inside(x, (&self.breaks_x[*idx], &self.breaks_x[*idx+1]));
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
				println!("curr = {:?}", &curr);
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
		println!("index = {:?}, x = {:?}, y0 = {:?}, y1 = {:?}", &index, &x, &y0, &y1);
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
		let c1 = coeffs[0]*coeffs[0]*xs;
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

impl PPData<f64>{
	
	pub fn intersection_with(&self, other: &PPData<f64>, index1: usize, index2: usize)->Option<(f64,f64)>{
		let mut coeffs = vec![0.0;4];
		coeffs[0] = self.coeffs[index1][0]-other.coeffs[index2][0];
		coeffs[1] = self.coeffs[index1][1]-other.coeffs[index2][1];
		coeffs[2] = self.coeffs[index1][2]-other.coeffs[index2][2];
		coeffs[3] = self.coeffs[index1][3]-other.coeffs[index2][3];
		println!("calculate_root = {:?}", &calculate_root(&coeffs));
		let x1 = self.breaks_x[index1] + calculate_root(&coeffs);
		let y1 = self.interpolate_for_index(&x1, index1)?;
		return Some((x1,y1));
	}
	
	pub fn intersection_with1(&self, other: &PPData<f64>, index1: usize, index2: usize)->Option<(f64,f64)>{
		println!("intersection_with1");
		let mut x1 = self.breaks_x[index1];
		let mut x2 = other.breaks_x[index2];
		let func = |x| self.interpolate(&x).unwrap()-other.interpolate(&x).unwrap();
		let mut y1 = func(x1);
		let mut y2 = func(x2);
		while (x1-x2).abs() > 1.0e-6 {
			println!("[{:?}-{:?}]", &x1, &x2);
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


/********************************************************************************************************************/

/// Multiple variables 1 degree of freedom interpolator. In some applications, instead of xx and yy pair, one might have multiple 1D variables xx, yy, zz, uu, vv, etc. If you select a principal variable tt, instead of constructing all possible pairs of variables, one can make only 2n-2 pairs (tt, xx), (xx, tt) to handle all possible cubic splines between the variables.
#[derive(Debug)]
pub struct MPPData<K,T>
where T: Float + fmt::Debug,
{
	pub keys: HashMap<K,Option<usize>>,  // variables are identified by keys, each keys corresponds to an index in pps vector, None for the principal variable
	pub pps: Vec<(PPData<T>,PPData<T>)>, // to? principal variable, from? principal variable Piecewise-Polynomials
	pub tt: Vec<T>,
}

impl<K,T: Float> MPPData<K,T>
where K : Eq + Hash, T : Float + std::fmt::Debug,
	{
	/// Initialize the interpolation structure using x0 variable as a principal variable
	pub fn new(key: K, tt: &[T])->Self{
		let mut keys: HashMap<K,Option<usize>> = HashMap::new();
		keys.insert(key, None); // initialize with the principal variable key + break point data
		let pps : Vec<(PPData<T>,PPData<T>)> = Vec::new();
		return Self{
			keys: keys,
			pps: pps,
			tt: tt.iter().map(|&idx| idx).collect(),
		};
	}
	/// Add another variable, key and values as a slice
	pub fn add_variable(&mut self, key: K, xx: &[T])->bool{
		match self.keys.get(&key){
			Some(Some(_)) => {
				return false; // cannot change an existing variable
			}
			Some(None) => {
				return false; // cannot change the principal variable
			}
			None => {
				// add a new variable
				let pp0 = makima(xx, &self.tt);
				let pp1 = makima(&self.tt, xx);
				let idx = self.pps.len();
				self.pps.push((pp0,pp1));
				self.keys.insert(key, Some(idx));
				return true;
			}
		}
	}
	
	/// interpolate between xi and yi (indicated by the keys) for an xi value
	pub fn interpolate<Q: ?Sized>(&self, keyx: &Q, keyy: &Q, x: &T)->Option<T>
	where K: Borrow<Q>,
		  Q : Hash + Eq,
	{
		match (self.keys.get(keyx), self.keys.get(keyy)){
			(None,_) | (_,None) => {return None;}
			(Some(Some(idx)), Some(None)) => {
				return self.pps[*idx].0.interpolate(x);
			}
			(Some(None), Some(Some(idx))) => {
				return self.pps[*idx].1.interpolate(x);
			}
			(Some(Some(idx0)), Some(Some(idx1)))=> {
				let t = self.pps[*idx0].0.interpolate(x)?;
				return  self.pps[*idx1].1.interpolate(&t);
			}
			(Some(None),Some(None)) => {return Some(*x);}
		}
	}
	
	pub fn interpolate_for_index<Q: ?Sized>(&self, keyx: &Q, keyy: &Q, x: &T, index: usize)->Option<T>
	where K: Borrow<Q>, Q: Hash + Eq,
	{
		match (self.keys.get(keyx), self.keys.get(keyy)){
			(None,_) | (_,None) => {return None;}
			(Some(Some(idx)), Some(None)) => {
				return self.pps[*idx].0.interpolate_for_index(x, index);
			}
			(Some(None), Some(Some(idx))) => {
				return self.pps[*idx].1.interpolate_for_index(x, index);
			}
			(Some(Some(idx0)), Some(Some(idx1)))=> {
				let t = self.pps[*idx0].0.interpolate_for_index(x, index)?;
				return  self.pps[*idx1].1.interpolate_for_index(&t, index);
			}
			(Some(None),Some(None)) => {return Some(*x);}
		}
	}
	
	pub fn interpolate_for_index_by_idx2pc(&self, idx: usize, x: &T, index: usize)->Option<T>{
		return self.pps[idx].0.interpolate_for_index(x,index);
	}
	
	pub fn interpolate_for_index_by_pc2idx(&self, idx: usize, x: &T, index: usize)->Option<T>{
		return self.pps[idx].1.interpolate_for_index(x,index);
	}
	
	pub fn get_break_for_index_by_key<Q: ?Sized>(&self, key: &Q, index: usize)->Option<T>
	where K: Borrow<Q>, Q: Hash + Eq,
	{
		match self.keys.get(key) {
			Some(Some(idx)) => {
				if index >= self.tt.len(){return None;}
				return Some(self.pps[*idx].0.breaks_x[index]);
			}
			Some(None) => {
				if index >= self.tt.len(){return None;}
				return Some(self.tt[index]);
			}
			None => {
				return None;
			}
		}
	}
	
	pub fn get_break_for_index_by_idx(&self, idx: usize, index: usize)->Option<T>{
		if index >= self.tt.len() && idx >= self.pps.len() {return None;}
		//println!("breaks = {:?}", self.pps[idx].0.breaks);
		return Some(self.pps[idx].0.breaks_x[index]);
	}
	
	pub fn number_variables(&self)->usize {
		return self.pps.len();
	}
	
	pub fn len(&self)->usize {
		return self.tt.len();
	}
	
	pub fn index(&self, key: K, x: &T)->Option<usize>{
		todo!();
	}
	
	pub fn index_by_idx(&self, idx: usize, x: &T)->Option<usize>{
		return self.pps[idx].0.index(x);
	}
	
}

/// Make interpolation structure from csv file contents
pub fn load_mpp_from_csv(datafile: &str)->Result<MPPData<String,f64>, Box<dyn Error>>{
	let file = File::open(&Path::new(datafile))?;
	let mut reader = Reader::from_reader(file);
	let headers = reader.headers()?.clone();
	let mut vecs : Vec<Vec<f64>> = (0..headers.len()).map(|_| Vec::new()).collect();
	for result in reader.records(){
		let record = result?;
		record.into_iter().enumerate().for_each(|(idx, strval)| vecs[idx].push(strval.parse::<f64>().unwrap_or(f64::NAN)));
	}
	
	let mut mpp : MPPData<String,f64> = MPPData::new(String::from(&headers[0]), &vecs[0]);
	
	for idx in 1..headers.len(){
		mpp.add_variable(String::from(&headers[idx]), &vecs[idx]);
	}
	
	return Ok(mpp);
}