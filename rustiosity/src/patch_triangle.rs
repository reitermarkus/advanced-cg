use crossbeam::epoch::{self, Atomic, Owned};

use std::sync::atomic::Ordering::{Acquire, Release};

use triangular::Triangular;
use vector::Vector;
use color::Color;

pub struct PatchTriangle {
  pub a: Vector,
  pub b: Vector,
  pub c: Vector,

  pub normal: Vector,
  pub color: Atomic<Color>,

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

    PatchTriangle { a: a, b: b, c: c, normal: normal, color: Atomic::new(Color::zero()), area: area }
  }

  pub fn set_color(&self, color: Color) {
    let mut color = Some(Owned::new(color));

    let guard = epoch::pin();

    loop {
      let old_color = self.color.load(Acquire, &guard);

      match self.color.cas(old_color, color, Release) {
        Ok(()) => return,
        Err(c) => color = c,
      }
    }
  }

  pub fn get_color(&self) -> Color {
    let guard = epoch::pin();

    loop {
      match self.color.load(Acquire, &guard) {
        Some(color) => return **color,
        None => continue,
      }
    }
  }
}

unsafe impl Send for PatchTriangle {}
unsafe impl Sync for PatchTriangle {}
