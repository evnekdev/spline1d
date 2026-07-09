// catmullrom.rs

use spline1d::{catmullrom, catmullrom_single_middle, catmullrom_single_middle_alpha};

fn eval(coeffs: [f64; 4], x_left: f64, x: f64) -> f64 {
    let dx = x - x_left;
    return coeffs.iter().fold(0.0, |acc, c| dx * acc + *c);
}

fn eval_alpha(y1: f64, y2: f64, alpha: [f64; 2], t: f64) -> f64 {
    return y1 * (1.0 - t) + y2 * t + (1.0 - t) * t * (alpha[0] + alpha[1] * t);
}

fn main() {
    let xx = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let yy = vec![0.0, 1.0, 1.5, 1.75, 2.0];

    let pp = catmullrom(&xx, &yy);

    for x in [0.0, 0.5, 1.0, 1.5, 2.5, 4.0] {
        println!("x = {x}, y = {:?}", pp.interpolate(&x));
    }

    let coeffs = catmullrom_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let alpha = catmullrom_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);

    println!("single middle at x=1.5: {}", eval(coeffs, 1.0, 1.5));
    println!("single middle alpha at t=0.5: {}", eval_alpha(1.0, 1.5, alpha, 0.5));
}
