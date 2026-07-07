// csvdata.rs
use std::mem::{size_of_val};
use spline1d::{binary_search_interval, PPData, makima, load_mpp_from_csv, SearchTree};

pub fn main(){
	let datafile = r"c:\_WORK\Code\Rust\workspace\spline1d\data\liq-mono.csv";
	let mpp = load_mpp_from_csv(datafile).unwrap();
	//println!("mpp = \n{:?}\n", &mpp);
	//println!("occupied size: {:?}", size_of_val(&mpp));
	let xmono = mpp.interpolate("T[C]", "xCaO(Liq)", &2500.0);
	let tliq  = mpp.interpolate("xCaO(Liq)", "T[C]", &0.95);
	println!("xmono = {:?}", &xmono);
	println!("tliq  = {:?}", &tliq);
}