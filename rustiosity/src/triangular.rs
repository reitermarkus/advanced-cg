extern crate rand;

use self::rand::distributions::{IndependentSample, Range};

use std::f64::EPSILON;

use ray::Ray;
use vector::Vector;

pub trait Triangular {
  fn a(&self) -> Vector;
  fn b(&self) -> Vector;
  fn c(&self) -> Vector;

  fn area(&self) -> f64;
  fn normal(&self) -> Vector;

  fn calculate_area(a: &Vector, b: &Vector, c: &Vector) -> f64 {
    ((b - a).cross_product(&(c - a)) / 2.0).length()
  }

  fn calculate_face_normal(a: &Vector, b: &Vector, c: &Vector) -> Vector {
    (b - a).cross_product(&(c - a)).normalize()
  }

  fn intersect(&self, ray: &Ray) -> f64 {
    let edge_1: Vector = &self.b() - &self.a();
    let edge_2: Vector = &self.c() - &self.a();

    let h: Vector = ray.dir.cross_product(&edge_2);
    let a: f64 = edge_1.dot_product(&h);

    if a > -EPSILON && a < EPSILON {
      return 0.0;
    }

    let f: f64 = 1.0 / a;
    let s: Vector = &ray.org - &self.a();
    let u: f64 = f * s.dot_product(&h);

    if u < 0.0 || u > 1.0 {
      return 0.0;
    }

    let q: Vector = s.cross_product(&edge_1);
    let v: f64 = f * ray.dir.dot_product(&q);

    if v < 0.0 || u + v > 1.0 {
      return 0.0;
    }

    let t: f64 = f * edge_2.dot_product(&q);

    if t <= EPSILON {
      return 0.0;
    }

    t
  }

  fn barycentric_coordinates_at(&self, p: &Vector) -> Vector {
    let area_p_b_c = ((&self.b() - p).cross_product(&(&self.c() - p)) / 2.0).length();
    let area_p_b_a = ((&self.b() - &self.a()).cross_product(&(p - &self.a())) / 2.0).length();

    let alpha = area_p_b_c / self.area();
    let beta  = area_p_b_a / self.area();
    let gamma = 1.0 - alpha - beta;

    Vector::new(alpha, beta, gamma)
  }

  fn random_sample(&self) -> Vector {
    let between = Range::new(0.0, 1.0);
    let mut rng = rand::thread_rng();

    let epsilon_0 = between.ind_sample(&mut rng) as f64;
    let epsilon_1 = between.ind_sample(&mut rng) as f64;

    let delta_0 = 1.0 - epsilon_0.sqrt();
    let delta_1 = epsilon_1 * epsilon_0.sqrt();
    let delta_2 = 1.0 - delta_0 - delta_1; // = sqrt(epsilon_0) * (1.0 - epsilon_1)

    let q: Vector = delta_0 * self.a() + delta_1 * self.b() + delta_2 * self.c();
    q
  }
}
