extern crate ordered_float;

use self::ordered_float::NotNaN;

use std::ops;
use std::fmt;

#[derive(Hash, PartialEq, Eq, Clone, Copy)]
pub struct Vector {
  pub x: NotNaN<f64>,
  pub y: NotNaN<f64>,
  pub z: NotNaN<f64>,
}

impl Vector {
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

impl fmt::Debug for Vector {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "Vector({}, {}, {})", self.x, self.y, self.z)
  }
}

impl ops::Add for Vector {
  type Output = Vector;

  fn add(self, other: Vector) -> Vector {
    Vector {
      x: self.x + other.x,
      y: self.y + other.y,
      z: self.z + other.z,
    }
  }
}

impl<'a> ops::Add<&'a Vector> for Vector {
  type Output = Vector;

  fn add(self, other: &'a Vector) -> Vector {
    Vector::new(self.x + other.x, self.y + other.y, self.z + other.z)
  }
}

impl ops::AddAssign for Vector {
  fn add_assign(&mut self, other: Vector) {
    *self = Vector::new(self.x + other.x, self.y + other.y, self.y + other.z);
  }
}

impl<'a, 'b> ops::Sub<&'b Vector> for &'a Vector {
  type Output = Vector;

  fn sub(self, other: &'b Vector) -> Vector {
    Vector::new(self.x - other.x, self.y - other.y, self.z - other.z)
  }
}

impl ops::Mul<f64> for Vector {
  type Output = Vector;
  fn mul(self, c: f64) -> Vector {
    Vector::new(self.x * c, self.y * c, self.z * c)
  }
}

impl ops::Mul<Vector> for i64 {
  type Output = Vector;
  fn mul(self, other: Vector) -> Vector {
    Vector::new(other.x * self as f64, other.y * self as f64, other.z * self as f64)
  }
}

impl ops::Mul<Vector> for f64 {
  type Output = Vector;
  fn mul(self, other: Vector) -> Vector {
    Vector::new(other.x * self, other.y * self, other.z * self)
  }
}

impl ops::Div<f64> for Vector {
  type Output = Vector;

  fn div(self, c: f64) -> Vector {
    Vector::new(self.x / c, self.y / c, self.z / c)
  }
}

impl<'a> ops::Div<f64> for &'a Vector {
  type Output = Vector;

  fn div(self, c: f64) -> Vector {
    Vector::new(self.x / c, self.y / c, self.z / c)
  }
}

impl ops::DivAssign<f64> for Vector {
  fn div_assign(&mut self, c: f64) {
    *self = Vector::new(self.x / c, self.y / c, self.z / c)
  }
}
