use spline1d::{
    makima_single_left, makima_single_middle, makima_single_middle_alpha, makima_single_right,
    pchip_single_left, pchip_single_middle, pchip_single_middle_alpha, pchip_single_right,
};

fn eval(coeffs: [f64; 4], x_left: f64, x: f64) -> f64 {
    let dx = x - x_left;
    return coeffs.iter().fold(0.0, |acc, c| dx * acc + *c);
}

fn main() {
    let c = makima_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = makima_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("makima middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("makima middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let c = pchip_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = pchip_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("pchip middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("pchip middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let _ = makima_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = makima_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = pchip_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = pchip_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
}

fn eval_alpha(y1: f64, y2: f64, alpha: [f64; 2], t: f64) -> f64 {
    return y1 * (1.0 - t) + y2 * t + (1.0 - t) * t * (alpha[0] + alpha[1] * t);
}
