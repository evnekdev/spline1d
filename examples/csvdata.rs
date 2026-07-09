// csvdata.rs
use spline1d::*;

pub fn main() {
    let datafile = r"c:\_WORK\Code\Rust\workspace\spline1d\data\liq-mono.csv";
    let mpp = load_multispline_from_csv(datafile).unwrap();
    //println!("mpp = \n{:?}\n", &mpp);
    //println!("occupied size: {:?}", size_of_val(&mpp));
    let xmono = mpp.interpolate("T[C]", "xCaO(Liq)", &2500.0);
    let tliq = mpp.interpolate("xCaO(Liq)", "T[C]", &0.95);
    println!("xmono = {:?}", &xmono);
    println!("tliq  = {:?}", &tliq);
}
