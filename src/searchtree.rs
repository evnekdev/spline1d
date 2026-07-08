// search_tree.rs

//! This module handles inverse cubic interpolation (finding all xs corresponding to a given value of y). Since in general case, there might be multiple values, the standard bincode which works for monotonous functions, will not work. Instread, a search tree is constructed which map monotonous regions of y->x mapping and root search happens in the tree subdivisions.
//! 
//! TODO - finish docs.

use std::fmt::{Debug};
use std::collections::{HashMap,BTreeSet};
use std::hash::{Hash};
use std::borrow::{Borrow};

use num::{Float};

use crate::multispline::MultiSpline;
use crate::binary_search_interval;

/*************************************************************************************************************/
/*************************************************************************************************************/

/// Search Node
#[derive(Debug, Clone)]
pub struct SearchNode<T>
where T: Float + Debug,
{
	interval: [usize;2],
	children: [usize;2],
	parent: usize,
	minmaxs: Vec<[T;2]>,
}

impl<T: Float + Debug> SearchNode<T>{
	
	pub fn root<K: Eq + Hash>(pps: &MultiSpline<K,T>) -> Self {
		let first : usize = 0;
		let last  : usize = pps.len()-1;
		let mut minmaxs : Vec<[T;2]> = Vec::new();
		for k in 0..pps.number_variables(){
			let first_val = pps.get_break_for_index_by_idx(k, first).unwrap();
			let last_val  = pps.get_break_for_index_by_idx(k, last).unwrap();
			minmaxs.push([T::min(first_val,last_val),T::max(first_val, last_val)]);
		}
		return Self {
			interval: [first,last],
			children: [0;2],
			parent: 0,
			minmaxs: minmaxs,
		};
	}
	
	pub fn null<K: Eq + Hash>(pps: &MultiSpline<K,T>) -> Self {
		let mut minmaxs : Vec<[T;2]> = Vec::new();
		for k in 0..pps.number_variables(){
			minmaxs.push([T::zero();2]);
		}
		return Self {
			interval: [0;2],
			children: [0;2],
			parent: 0,
			minmaxs: minmaxs,
		};
	}
	
}

/*************************************************************************************************************/
/*************************************************************************************************************/

#[derive(Debug)]
pub struct SearchTree<'a, K, T>
where T: Float + Debug, K : Eq + Hash,
{
	pub nodes: Vec<SearchNode<T>>,
	pub pps: &'a MultiSpline<K,T>,
}

impl<'a, K,T> SearchTree<'a, K,T>
where T: Float + Debug, K : Eq + Hash,
{
	pub fn new(pps: &'a MultiSpline<K,T>)->Self {
		let mut tree = Self {
			nodes: vec![SearchNode::null(pps), SearchNode::root(pps)],
			pps: pps,
		};
		let extrema = tree.search_extrema_linear();
		//println!("extrema = {:?}", &extrema);
		for k in 0..extrema.len(){
			tree.split_node_at(1, extrema[k].0);
		}
		return tree;
	}
	
	pub fn search_extrema_linear(&self)->Vec<(usize,usize)>{ // index, variable #
		let mut previous : Vec<T> = vec![T::zero();self.nodes[1].minmaxs.len()];
		let mut indices : Vec<usize> = vec![0usize;self.nodes[1].minmaxs.len()];
		let mut setprev : Vec<bool>  = vec![false;self.nodes[1].minmaxs.len()];
		let mut increasing : Vec<bool> = vec![false;self.nodes[1].minmaxs.len()];
		let mut res : Vec<(usize,usize)> = Vec::new();
		for k in 0..self.pps.len(){
			for m in 0..indices.len(){
				let brk = self.pps.get_break_for_index_by_idx(m, k).unwrap();
				if k == 0 {
					previous[m] = brk;
					continue;
				}
				if brk < previous[m] {
					if k >= 2 && increasing[m] {
						res.push((k-1,m));
					}
					indices[m] = k;
					previous[m] = brk;
					increasing[m] = false;
					setprev[m] = true;
				} else {
					if k >= 2 && !increasing[m]{
						res.push((k-1,m));
					}
					indices[m] = k;
					previous[m] = brk;
					increasing[m] = true;
					setprev[m] = true;
				}
			}
		}
		return res;
	}

	pub fn split_node_at(&mut self, nindex: usize, new_index: usize)->Option<(usize,usize)>{
		if new_index <= self.nodes[nindex].interval[0] || new_index >= self.nodes[nindex].interval[1] {return None;}
		if self.nodes[nindex].children[0] > 0 {
			match self.split_node_at(self.nodes[nindex].children[0], new_index){
				Some(res) => return Some(res),
				None => {},
			}
			match self.split_node_at(self.nodes[nindex].children[1], new_index){
				Some(res) => return Some(res),
				None => return None,
			}
		}
		let mut left  = self.nodes[nindex].clone();
		let mut right = self.nodes[nindex].clone();
		let idx_left = self.nodes.len();
		let idx_right = self.nodes.len()+1;
		self.nodes.push(left);
		self.nodes.push(right);
		self.nodes[nindex].children = [idx_left,idx_right];
		self.nodes[idx_left].parent = nindex;
		self.nodes[idx_right].parent = nindex;
		for m in 0..self.nodes[nindex].minmaxs.len(){
			let interval_left  = self.nodes[nindex].interval[0];
			let interval_right = self.nodes[nindex].interval[1];
			let value_left  = self.pps.get_break_for_index_by_idx(m, interval_left).unwrap();
			let value_right = self.pps.get_break_for_index_by_idx(m, interval_right).unwrap();
			let value_middle = self.pps.get_break_for_index_by_idx(m, new_index).unwrap();
			self.nodes[idx_left].minmaxs[m]  = [T::min(value_left, value_middle), T::max(value_left, value_middle)];
			self.nodes[idx_right].minmaxs[m] = [T::min(value_middle, value_right),T::max(value_middle, value_right)];
		}
		
		self.nodes[idx_left].interval  = [self.nodes[nindex].interval[0], new_index];
		self.nodes[idx_right].interval = [new_index, self.nodes[nindex].interval[1]];
		self.update_parent(idx_left);
		self.update_parent(idx_right);
		return Some((idx_left,idx_right));
	}
	
	fn update_parent(&mut self, nindex: usize){
		if nindex == 0 {return;}
		let parent = self.nodes[nindex].parent;
		if parent == 0 {return;}
		for k in 0..self.nodes[nindex].minmaxs.len(){
			// check minimum value
			if self.nodes[nindex].minmaxs[k][0] < self.nodes[parent].minmaxs[k][0]{
				self.nodes[parent].minmaxs[k][0] = self.nodes[nindex].minmaxs[k][0];
			}
			// check maximum value
			if self.nodes[nindex].minmaxs[k][1] > self.nodes[parent].minmaxs[k][1]{
				self.nodes[parent].minmaxs[k][1] = self.nodes[nindex].minmaxs[k][1];
			}
		}
	}
	
	/// TODO
	pub fn interval_indices<Q: ?Sized>(&self, key: &Q, x: &T)->Vec<T>
	where K: Borrow<Q>, Q : Hash + Eq,
	{
		match self.pps.keys.get(key){
			None => return Vec::new(),
			Some(None) => {
				todo!();
			}
			Some(Some(idx)) => {
				todo!();
			}
		}
	}
	
	pub fn interval_indices_by_idx(&self, idx: usize, x: &T)->BTreeSet<usize>{
		let mut res : BTreeSet<usize> = BTreeSet::new();
		self.interval_indices_by_idx_(idx, x, &mut res, 1);
		return res;
	}
	
	fn interval_indices_by_idx_(&self, idx: usize, x: &T, res: &mut BTreeSet<usize>, index: usize){
		if *x < self.nodes[index].minmaxs[idx][0] || *x > self.nodes[index].minmaxs[idx][1] {
			return;
		}
		if self.nodes[index].children[0] > 0 {
			self.interval_indices_by_idx_(idx, x, res, self.nodes[index].children[0]);
			self.interval_indices_by_idx_(idx, x, res, self.nodes[index].children[1]);
			return;
		}
		match self.binary_search_monotonic_by_idx(idx, x, self.nodes[index].interval[0], self.nodes[index].interval[1]){
			Some(ix) => {
				res.insert(ix);
			}
			None => {}
		}
	}
	
	fn binary_search_monotonic_by_idx(&self, idx: usize, x: &T, index0: usize, index1: usize)->Option<usize>{
		return binary_search_interval(index1 - index0 + 1, x, |ix| self.pps.pps[idx].0.breaks_x[ix + index0]).map(|idx| idx + index0);
	}
	
	fn binary_search_monotonic_principal(&self, x: &T, index0: usize, index1: usize)->Option<usize>{
		return binary_search_interval(index1 - index0 + 1, x, |ix| self.pps.tt[ix + index0]).map(|idx| idx + index0);
	}
	
	pub fn interpolate<Q: ?Sized>(&self, keyx: &Q, keyy: &Q, x: &T)->Vec<T>
	where K: Borrow<Q>, Q : Hash + Eq,
	{
		match (self.pps.keys.get(keyx), self.pps.keys.get(keyy)){
			(None, _) | (_, None) => {return Vec::new();}
			(Some(None), Some(Some(idx))) => {
				match self.binary_search_monotonic_principal(x, self.nodes[1].interval[0], self.nodes[1].interval[1]){
					Some(index) => {
						return vec![self.pps.interpolate_for_index_by_pc2idx(*idx, x, index).unwrap()];
					}
					None => {
						return vec![];
					}
				}
				return self.pps.pps[*idx].0.interpolate(x).into_iter().collect();
			}
			(Some(Some(idx)), Some(None)) => {
				let indices = self.interval_indices_by_idx(*idx, x);
				let mut res : Vec<T> = Vec::new();
				for index in indices.into_iter(){
					res.push(self.pps.interpolate_for_index_by_idx2pc(*idx, x, index).unwrap());
				}
				return res;
			}
			(Some(Some(idx0)), Some(Some(idx1))) => {
				let indices = self.interval_indices_by_idx(*idx0, x);
				let mut res : Vec<T> = Vec::new();
				for index in indices.into_iter(){
					let t   = self.pps.interpolate_for_index_by_idx2pc(*idx0, x, index).unwrap();
					res.push(self.pps.interpolate_for_index_by_pc2idx(*idx1, &t, index).unwrap());
				}
				return res;
			}
			(Some(None), Some(None)) => {return vec![*x];}
		}
	}
	
}

/*************************************************************************************************************/
/*************************************************************************************************************/