// main.rs (splines library)

use std::mem::{size_of_val};

use spline1d::{binary_search_interval, PPData, makima, load_mpp_from_csv, SearchTree};

pub fn main(){
	/*
	let vec : Vec<f64> = vec![0.0, 0.1, 0.25, 0.32, 0.5];
	let diffs1 : Vec<f64> = diff(&vec).collect();
	let diffs2 : Vec<f64> = kernel_conv(&vec, &[-1.0, 1.0]).collect();
	let diffdiff: Vec<f64> = kernel_conv(&vec, &[-1.0, 2.0, -1.0]).collect();
	//println!("diffs1 = {:?}", &diffs1);
	//println!("diffs2 = {:?}", &diffs2);
	//println!("diffdiff = {:?}", &diffdiff);
	let sval = 0.0;
	let interval = binary_search_interval(vec.len(), &sval, |loc| slice_locator(&vec, loc));
	//println!("interval = {:?}", &interval);
	let pp = PPData{breaks: vec![0.0,0.2,0.4,0.55,0.65],
					coeffs: vec![[1074.93702171841,-1247.84215434368,-579.278039999999,2845.16205000000],
								 [-3310.42265949697,-585.757622444289,-949.422459131263,2687.99225200000],
								 [11715.1755049226,-6054.47449963650,-1580.97622724862,2448.19407400000],
								 [12313.6022776045,-5823.54912218758,-2606.54423055729,2114.36068100000]]};
	let xx : Vec<f64> = vec![0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55];
	for x in xx.iter(){
		let y = pp.interpolate(x);
		println!("{:?},{:?}", x, &y);
	}
	*/
	test_makima();
	//test_mpp_csv();
	//test_search_tree();
}


pub fn test_makima(){
	let xx = vec![0.0, 0.2, 0.4, 0.55, 0.65];
	let yy = vec![2845.2, 2688.0, 2448.2, 2114.4, 1807.8];
	let pp = makima(&xx, &yy);
	let xx : Vec<f64> = vec![0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55];
	for x in xx.iter(){
		let y = pp.interpolate(x);
		println!("{:?},{:?}", x, &y);
	}
}
/*
1074.93702171841	-1247.84215434368	-579.278039999999	2845.16205000000
-3310.42265949697	-585.757622444289	-949.422459131263	2687.99225200000
11715.1755049226	-6054.47449963650	-1580.97622724862	2448.19407400000
12313.6022776045	-5823.54912218758	-2606.54423055729	2114.36068100000
*/

pub fn test_mpp_csv(){
	let datafile = r"c:\_WORK\Code\Rust\workspace\splines\data\liq-mono.csv";
	let mpp = load_mpp_from_csv(datafile).unwrap();
	//println!("mpp = \n{:?}\n", &mpp);
	//println!("occupied size: {:?}", size_of_val(&mpp));
	let xmono = mpp.interpolate("T[C]", "xCaO(Liq)", &2500.0);
	let tliq  = mpp.interpolate("xCaO(Liq)", "T[C]", &0.95);
	println!("xmono = {:?}", &xmono);
	println!("tliq  = {:?}", &tliq);
}

pub fn test_search_tree(){
	let datafile = r"c:\_WORK\Code\Rust\workspace\splines\data\poly.csv";
	let mpp = load_mpp_from_csv(datafile).unwrap();
	let mut tree = SearchTree::new(&mpp);
	let extrema = tree.search_extrema_linear();
	for k in 0..extrema.len(){
		tree.split_node_at(1, extrema[k].0);
	}
	let x1 = 0.10;
	let indices = tree.interval_indices_by_idx(0,&x1);
	let values  = tree.interpolate("x1", "x2", &x1);
	println!("{:?}", &tree.nodes);
	println!("{:?}", &tree.pps);
	println!("{:?}", &extrema);
	println!("indices = {:?}", &indices);
	println!("values  = {:?}", &values);
}