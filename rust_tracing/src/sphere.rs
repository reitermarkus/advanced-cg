use std::any::Any;

use scene_object::SceneObject;
use scene_object::ReflType;
use vector::Vector;
use color::Color;
use ray::Ray;

pub struct Sphere {
  pub radius: f64,
  pub position: Vector,

  pub emission: Vector,
  pub color: Vector,

  pub refl: ReflType
}

impl SceneObject for Sphere {
  fn color(&self) -> Color { self.color }
  fn emission(&self) -> Color { self.emission }
  fn refl(&self) -> ReflType { self.refl }
  fn is_sphere(&self) -> bool { true }

  fn as_any(&self) -> &dyn Any {
    self
  }

  fn intersect(&self, ray: &Ray) -> f64 {
    /* Check for ray-sphere intersection by solving for t:
      t^2*d.d + 2*t*(o-p).d + (o-p).(o-p) - R^2 = 0 */
    let op = self.position - ray.org;
    let eps = 1e-4;
    let b = op.dot_product(&ray.dir);
    let mut radicant = b * b - op.dot_product(&op) + self.radius * self.radius;

    /* No intersection */
    if radicant < 0.0 {
      return 0.0;
    }

    radicant = radicant.sqrt();

    /* Check smaller root first */
    let mut t = b - radicant;

    if t > eps {
      return t;
    }

    t = b + radicant;

    /* Check second root */
    if t > eps  {
      return t;
    }

    /* No intersection in ray direction */
    return 0.0;
  }
}

impl Sphere {
  pub fn new(radius: f64, position: Vector, emission: Vector, color: Vector, refl: ReflType) -> Sphere {
    Sphere {
      radius: radius,
      position: position,
      emission: emission,
      color: color,
      refl: refl
    }
  }
}
