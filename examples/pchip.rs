// pchip.rs

use spline1d::pchip;

fn main() {
    let xx = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let yy = vec![0.0, 1.0, 1.5, 1.75, 2.0];

    let pp = pchip(&xx, &yy);

    for x in [0.0, 0.5, 1.0, 1.5, 2.5, 4.0] {
        println!("x = {x}, y = {:?}", pp.interpolate(&x));
    }
}
