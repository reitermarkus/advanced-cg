extern crate ordered_float;

use self::ordered_float::NotNaN;

use std::ops::{Add, AddAssign, Sub, Mul, Div, DivAssign};
use std::fmt::{Debug, Formatter, Result};
use std::iter::Sum;

#[derive(Hash, PartialEq, Eq, Clone, Copy)]
pub struct Vector {
  pub x: NotNaN<f64>,
  pub y: NotNaN<f64>,
  pub z: NotNaN<f64>,
}

impl Vector {
  pub fn zero() -> Vector {
    Vector::new(0.0, 0.0, 0.0)
  }

  pub fn new<T: Into<NotNaN<f64>>>(x: T, y: T, z: T) -> Vector {
    Vector { x: x.into(), y: y.into(), z: z.into() }
  }

  pub fn entrywise_product(&self, other: &Vector) -> Vector {
    Vector {
      x: self.x * other.x,
      y: self.y * other.y,
      z: self.z * other.z,
    }
  }

  pub fn dot_product(&self, other: &Vector) -> f64 {
    (self.x * other.x + self.y * other.y + self.z * other.z).into_inner()
  }

  pub fn cross_product(&self, other: &Vector) -> Vector {
    Vector {
      x: (self.y * other.z) - (self.z * other.y),
      y: (self.z * other.x) - (self.x * other.z),
      z: (self.x * other.y) - (self.y * other.x),
    }
  }

  pub fn distance(&self, other: &'static Vector) -> f64 {
    (self - other).length()
  }

  pub fn normalize(&self) -> Vector {
    self / self.length()
  }

  pub fn length_squared(&self) -> f64 {
    self.dot_product(self)
  }

  pub fn length(&self) -> f64 {
    self.length_squared().sqrt()
  }
}

impl Debug for Vector {
  fn fmt(&self, f: &mut Formatter) -> Result {
    write!(f, "Vector({}, {}, {})", self.x, self.y, self.z)
  }
}

impl Add for Vector {
  type Output = Vector;

  fn add(self, other: Vector) -> Vector {
    Vector {
      x: self.x + other.x,
      y: self.y + other.y,
      z: self.z + other.z,
    }
  }
}

impl<'a> Add<&'a Vector> for Vector {
  type Output = Vector;

  fn add(self, other: &'a Vector) -> Vector {
    Vector::new(self.x + other.x, self.y + other.y, self.z + other.z)
  }
}

impl AddAssign for Vector {
  fn add_assign(&mut self, other: Vector) {
    *self = *self + &other;
  }
}

impl<'a, 'b> Sub<&'b Vector> for &'a Vector {
  type Output = Vector;

  fn sub(self, other: &'b Vector) -> Vector {
    Vector::new(self.x - other.x, self.y - other.y, self.z - other.z)
  }
}

impl Mul<f64> for Vector {
  type Output = Vector;
  fn mul(self, c: f64) -> Vector {
    Vector::new(self.x * c, self.y * c, self.z * c)
  }
}

impl Mul<Vector> for i64 {
  type Output = Vector;
  fn mul(self, other: Vector) -> Vector {
    other * (self as f64)
  }
}

impl Mul<Vector> for f64 {
  type Output = Vector;
  fn mul(self, other: Vector) -> Vector {
    other * self
  }
}

impl Div<f64> for Vector {
  type Output = Vector;

  fn div(self, c: f64) -> Vector {
    Vector::new(self.x / c, self.y / c, self.z / c)
  }
}

impl<'a> Div<f64> for &'a Vector {
  type Output = Vector;

  fn div(self, c: f64) -> Vector {
    Vector::new(self.x / c, self.y / c, self.z / c)
  }
}

impl DivAssign<f64> for Vector {
  fn div_assign(&mut self, c: f64) {
    *self = Vector::new(self.x / c, self.y / c, self.z / c)
  }
}

impl Sum for Vector {
  fn sum<I>(iter: I) -> Vector where I: Iterator<Item = Self> {
    let mut zero = Vector::new(0.0, 0.0, 0.0);

    for vector in iter {
      zero += vector;
    }

    zero
  }
}
