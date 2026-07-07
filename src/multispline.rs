// spline1d::multispline.rs


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
use crate::pchip::{pchip};
use crate::spline::Spline;
use crate::solve::{calculate_root};

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/

/// Multiple variables 1 degree of freedom interpolator. In some applications, instead of xx and yy pair, one might have multiple 1D variables xx, yy, zz, uu, vv, etc. If you select a principal variable tt, instead of constructing all possible pairs of variables, one can make only 2n-2 pairs (tt, xx), (xx, tt) to handle all possible cubic splines between the variables.
#[derive(Debug)]
pub struct MultiSpline<K,T>
where T: Float + fmt::Debug,
{
	pub keys: HashMap<K,Option<usize>>,  // variables are identified by keys, each keys corresponds to an index in pps vector, None for the principal variable
	pub pps: Vec<(Spline<T>,Spline<T>)>, // to? principal variable, from? principal variable Piecewise-Polynomials
	pub tt: Vec<T>,
}

impl<K,T: Float> MultiSpline<K,T>
where K : Eq + Hash, T : Float + std::fmt::Debug,
	{
	/// Initialize the interpolation structure using x0 variable as a principal variable
	pub fn new(key: K, tt: &[T])->Self{
		let mut keys: HashMap<K,Option<usize>> = HashMap::new();
		keys.insert(key, None); // initialize with the principal variable key + break point data
		let pps : Vec<(Spline<T>,Spline<T>)> = Vec::new();
		return Self {
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

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/

/// Make interpolation structure from csv file contents
pub fn load_mpp_from_csv(datafile: &str)->Result<MultiSpline<String,f64>, Box<dyn Error>>{
	let file = File::open(&Path::new(datafile))?;
	let mut reader = Reader::from_reader(file);
	let headers = reader.headers()?.clone();
	let mut vecs : Vec<Vec<f64>> = (0..headers.len()).map(|_| Vec::new()).collect();
	for result in reader.records(){
		let record = result?;
		record.into_iter().enumerate().for_each(|(idx, strval)| vecs[idx].push(strval.parse::<f64>().unwrap_or(f64::NAN)));
	}
	
	let mut mpp : MultiSpline<String,f64> = MultiSpline::new(String::from(&headers[0]), &vecs[0]);
	
	for idx in 1..headers.len(){
		mpp.add_variable(String::from(&headers[idx]), &vecs[idx]);
	}
	
	return Ok(mpp);
}

/*****************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************/