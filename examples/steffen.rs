use spline1d::{steffen, steffen_single_middle, steffen_single_middle_alpha};

fn main() {
    let xx = [0.0, 1.0, 2.0, 4.0];
    let yy = [0.0, 1.0, 1.5, 3.0];

    let spline = steffen(&xx, &yy);
    println!("y(1.5) = {:?}", spline.interpolate(&1.5));

    let coeffs = steffen_single_middle(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);
    let alpha = steffen_single_middle_alpha(xx[0], yy[0], xx[1], yy[1], xx[2], yy[2], xx[3], yy[3]);

    println!("local cubic coeffs = {:?}", coeffs);
    println!("alpha coeffs = {:?}", alpha);
}
