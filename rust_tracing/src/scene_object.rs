use color::Color;
use ray::Ray;

use std::any::Any;

#[derive(Copy, Clone)]
pub enum ReflType {
  DIFF,
  SPEC,
  REFR,
  GLOS,
  TRAN
}

pub trait SceneObject: Sync + Send {
  fn color(&self) -> Color;
  fn emission(&self) -> Color;
  fn refl(&self) -> ReflType;

  fn as_any(&self) -> &dyn Any;

  fn intersect(&self, ray: &Ray) -> f64;
}
