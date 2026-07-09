use spline1d::{
    akima_single_left, akima_single_middle, akima_single_middle_alpha, akima_single_right,
    cardinal_single_left, cardinal_single_middle, cardinal_single_middle_alpha, cardinal_single_right,
    catmullrom_single_left, catmullrom_single_middle, catmullrom_single_middle_alpha, catmullrom_single_right,
    fritschbutland_single_left, fritschbutland_single_middle, fritschbutland_single_middle_alpha, fritschbutland_single_right,
    makima_single_left, makima_single_middle, makima_single_middle_alpha, makima_single_right,
    pchip_single_left, pchip_single_middle, pchip_single_middle_alpha, pchip_single_right,
    steffen_single_left, steffen_single_middle, steffen_single_middle_alpha, steffen_single_right,
};

fn eval(coeffs: [f64; 4], x_left: f64, x: f64) -> f64 {
    let dx = x - x_left;
    return coeffs.iter().fold(0.0, |acc, c| dx * acc + *c);
}

fn eval_alpha(y1: f64, y2: f64, alpha: [f64; 2], t: f64) -> f64 {
    return y1 * (1.0 - t) + y2 * t + (1.0 - t) * t * (alpha[0] + alpha[1] * t);
}

fn main() {
    let c = akima_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = akima_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("akima middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("akima middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let c = makima_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = makima_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("makima middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("makima middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let c = pchip_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = pchip_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("pchip middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("pchip middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let c = steffen_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = steffen_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("steffen middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("steffen middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let c = catmullrom_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = catmullrom_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("catmullrom middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("catmullrom middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let tension = 0.25;
    let c = cardinal_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75, tension);
    let a = cardinal_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75, tension);
    println!("cardinal middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("cardinal middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let c = fritschbutland_single_middle(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    let a = fritschbutland_single_middle_alpha(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, 3.0, 1.75);
    println!("fritschbutland middle at 1.5 = {}", eval(c, 1.0, 1.5));
    println!("fritschbutland middle alpha at t=0.5 = {}", eval_alpha(1.0, 1.5, a, 0.5));

    let _ = akima_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = akima_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = makima_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = makima_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = pchip_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = pchip_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = steffen_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = steffen_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = catmullrom_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = catmullrom_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = cardinal_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, tension);
    let _ = cardinal_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5, tension);
    let _ = fritschbutland_single_left(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
    let _ = fritschbutland_single_right(0.0, 0.0, 1.0, 1.0, 2.0, 1.5);
}
