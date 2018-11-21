use std::f64::consts::PI;
use vector::Vector;
use rand::Rng;

pub fn drand48(from: f64, to: f64) -> f64 {
  rand::thread_rng().gen_range(from, to)
}
pub fn random_direction(direction: Vector, max_angle: f64)  -> Vector {
    /* Set up local orthogonal coordinate system u,v,w on surface */
  let w = direction;
  let u = if w.x.abs() > 0.1 { Vector::new(0.0, 1.0, 0.0) } else { Vector::new(1.0, 0.0, 0.0) };
  let u = u.cross_product(&w).normalize();

  let v = w.cross_product(&u);

  let eps1 = drand48(0.0, 1.0);
  let eps2 = drand48(0.0, 1.0);

  let cos_a = 1.0 - eps1 + eps1 * max_angle.cos();
  let sin_a = (1.0 - cos_a * cos_a).sqrt();
  let phi = 2.0 * PI * eps2;

  /* Random reflection vector d */
  let l = u * phi.cos() * sin_a +
          v * phi.sin() * sin_a +
          w * cos_a;

  l.normalize()
}