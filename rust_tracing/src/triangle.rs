use std::any::Any;

use color::Color;
use vector::Vector;
use scene_object::ReflType;
use scene_object::SceneObject;
use ray::Ray;

use std::f64::EPSILON;

pub struct Triangle {
  pub a: Vector,
  pub b: Vector,
  pub c: Vector,

  pub b_rel: Vector,
  pub c_rel: Vector,

  pub normal: Vector,

  pub emission: Vector,
  pub color: Vector,

  pub refl: ReflType
}

impl SceneObject for Triangle {
  fn color(&self) -> Color { self.color }
  fn emission(&self) -> Color { self.emission }
  fn refl(&self) -> ReflType { self.refl }
  fn is_sphere(&self) -> bool { false }

  fn as_any(&self) -> &dyn Any { self }

  fn intersect(&self, ray: &Ray) -> f64 {
    let edge_1: Vector = &self.b - &self.a;
    let edge_2: Vector = &self.c - &self.a;

    let h: Vector = ray.dir.cross_product(&edge_2);
    let a: f64 = edge_1.dot_product(&h);

    if a > -EPSILON && a < EPSILON {
      return 0.0;
    }

    let f: f64 = 1.0 / a;
    let s: Vector = &ray.org - &self.a;
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
}

impl Triangle {
  pub fn new(a: Vector, b_rel: Vector, c_rel: Vector, emission: Color, color: Color, refl: ReflType) -> Triangle {
    let b = a + b_rel;
    let c = a + c_rel;

    let normal = Self::calculate_face_normal(&a, &b, &c);

    Triangle {
      a: a, b: b, c: c,
      b_rel: b_rel, c_rel: c_rel,
      normal: normal,
      emission: emission, color: color,
      refl: refl
    }
  }

  fn calculate_face_normal(a: &Vector, b: &Vector, c: &Vector) -> Vector {
    (b - a).cross_product(&(c - a)).normalize()
  }
}
