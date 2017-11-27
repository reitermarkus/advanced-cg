use std::vec::Vec;

use color::Color;
use vector::Vector;
use triangular::Triangular;
use patch_triangle::PatchTriangle;

pub struct Triangle {
  pub a: Vector,
  pub b: Vector,
  pub c: Vector,

  pub b_rel: Vector,
  pub c_rel: Vector,

  pub normal: Vector,

  pub area: f64,

  pub emission: Vector,
  pub color: Vector,

  pub divisions: u64,

  pub patches: Vec<Color>,
  pub sub_triangles: Vec<PatchTriangle>,
}

impl Triangular for Triangle {
  fn a(&self) -> Vector { self.a }
  fn b(&self) -> Vector { self.b }
  fn c(&self) -> Vector { self.c }

  fn area(&self) -> f64 { self.area }
  fn normal(&self) -> Vector { self.normal }
}

impl Triangle {
  pub fn new(a: Vector, b_rel: Vector, c_rel: Vector, emission: Color, color: Color) -> Triangle {
    let b = a + b_rel;
    let c = a + c_rel;

    let area = Self::calculate_area(&a, &b, &c);
    let normal = Self::calculate_face_normal(&a, &b, &c);

    Triangle {
      a: a, b: b, c: c,
      b_rel: b_rel, c_rel: c_rel,
      normal: normal,
      area: area,
      emission: emission, color: color,
      divisions: 0, patches: Vec::new(), sub_triangles: Vec::new()
    }
  }

  pub fn init_patches(&mut self, divisions: u64) {
    self.divisions = divisions;
    self.patches = vec![Color::zero(); divisions.pow(2) as usize];
    self.divide(divisions)
  }

  fn divide(&mut self, divisions: u64) {
    let delta_x = self.b_rel / (divisions as f64);
    let delta_y = self.c_rel / (divisions as f64);

    self.sub_triangles = Vec::new();

    // Loop through patches, from bottom to top, left to right.
    let mut row_size = divisions as i64;
    let mut row = 0;
    while row_size > 0 {
      for x in 0..row_size {
        let y: i64 = row / 2;

        let offset = row % 2;

        let v1 = self.a + (x + offset) * delta_x + y * delta_y;
        let v2 = self.a + (x + 1) * delta_x + (y + offset) * delta_y;
        let v3 = self.a + x * delta_x + (y + 1) * delta_y;

        self.sub_triangles.push(PatchTriangle::new(v1, v2, v3))
      }

      row += 1;

      if row % 2 == 1 {
        row_size -= 1;
      }
    }
  }
}
