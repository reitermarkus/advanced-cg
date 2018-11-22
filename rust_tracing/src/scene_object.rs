use downcast_rs::Downcast;

use color::Color;
use ray::Ray;

#[derive(Copy, Clone)]
pub enum ReflType {
  DIFF,
  SPEC,
  REFR,
  GLOS,
  TRAN
}

pub trait SceneObject: Downcast + Sync + Send {
  fn color(&self) -> Color;
  fn emission(&self) -> Color;
  fn refl(&self) -> ReflType;

  fn intersect(&self, ray: &Ray) -> f64;
}

impl_downcast!(SceneObject);
