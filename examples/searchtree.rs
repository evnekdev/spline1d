// searchtree.rs
use spline1d::*;

pub fn main() {
    let datafile = r"c:\_WORK\Code\Rust\workspace\spline1d\data\poly.csv";
    let mpp = load_multispline_from_csv(datafile).unwrap();
    let mut tree = SearchTree::new(&mpp);
    let extrema = tree.search_extrema_linear();
    for k in 0..extrema.len() {
        tree.split_node_at(1, extrema[k].0);
    }
    let x1 = 0.10;
    let indices = tree.interval_indices_by_idx(0, &x1);
    let values = tree.interpolate("x1", "x2", &x1);
    println!("{:?}", &tree.nodes);
    println!("{:?}", &tree.pps);
    println!("{:?}", &extrema);
    println!("indices = {:?}", &indices);
    println!("values  = {:?}", &values);
}
