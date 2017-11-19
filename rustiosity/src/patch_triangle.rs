use triangular::Triangular;
use vector::Vector;

pub struct PatchTriangle {
  pub a: Vector,
  pub b: Vector,
  pub c: Vector,

  pub normal: Vector,

  pub area: f64,
}

impl Triangular for PatchTriangle {
  fn a(&self) -> Vector { self.a }
  fn b(&self) -> Vector { self.b }
  fn c(&self) -> Vector { self.c }

  fn area(&self) -> f64 { self.area }
  fn normal(&self) -> Vector { self.normal }
}

impl PatchTriangle {
  pub fn new(a: Vector, b: Vector, c: Vector) -> PatchTriangle {
    let area = Self::calculate_area(&a, &b, &c);
    let normal = Self::calculate_face_normal(&a, &b, &c);

    PatchTriangle { a: a, b: b, c: c, normal: normal, area: area }
  }
}
