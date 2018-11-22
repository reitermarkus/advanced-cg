use std::fmt::{Debug, Formatter, Result};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Sub, Mul, Div, DivAssign, Neg};

#[derive(PartialEq, Clone, Copy)]
pub struct Vector {
  pub x: f64,
  pub y: f64,
  pub z: f64,
}

impl Vector {
  pub fn zero() -> Self {
    Self::new(0.0, 0.0, 0.0)
  }

  pub fn new(x: f64, y: f64, z: f64) -> Self {
    Self { x: x, y: y, z: z }
  }

  pub fn entrywise_product(&self, other: Self) -> Self {
    Self {
      x: self.x * other.x,
      y: self.y * other.y,
      z: self.z * other.z,
    }
  }

  pub fn dot_product(&self, other: &Self) -> f64 {
    self.x * other.x + self.y * other.y + self.z * other.z
  }

  pub fn cross_product(&self, other: &Self) -> Self {
    Self {
      x: (self.y * other.z) - (self.z * other.y),
      y: (self.z * other.x) - (self.x * other.z),
      z: (self.x * other.y) - (self.y * other.x),
    }
  }

  pub fn distance(&self, other: &'static Self) -> f64 {
    (self - other).length()
  }

  pub fn normalize(&self) -> Self {
    self / self.length()
  }

  pub fn length_squared(&self) -> f64 {
    self.dot_product(self)
  }

  pub fn length(&self) -> f64 {
    self.length_squared().sqrt()
  }

  pub fn clamp(&self, min: f64, max: f64) -> Self {
    let clamp_sub = |input: f64, min: f64, max: f64| match input {
      i if i > max => max,
      i if i < min => min,
      i => i,
    };

    Self::new(
      clamp_sub(self.x, min, max),
      clamp_sub(self.y, min, max),
      clamp_sub(self.z, min, max),
    )
  }

  pub fn max(&self) -> f64 {
    let maxf = |a: f64, b: f64| if a > b { a } else { b };
    maxf(self.x, maxf(self.y, self.z))
  }
}

impl Debug for Vector {
  fn fmt(&self, f: &mut Formatter) -> Result {
    write!(f, "Vector({}, {}, {})", self.x, self.y, self.z)
  }
}

impl Add for Vector {
  type Output = Self;

  fn add(self, other: Self) -> Self::Output {
    Self::Output {
      x: self.x + other.x,
      y: self.y + other.y,
      z: self.z + other.z,
    }
  }
}

impl Sub for Vector {
  type Output = Self;

  fn sub(self, other: Self) -> Self::Output {
    Self::Output {
      x: self.x - other.x,
      y: self.y - other.y,
      z: self.z - other.z,
    }
  }
}

impl<'a> Add<&'a Vector> for Vector {
  type Output = Self;

  fn add(self, other: &'a Self) -> Self::Output {
    Self::Output::new(self.x + other.x, self.y + other.y, self.z + other.z)
  }
}

impl Add<f64> for Vector {
  type Output = Self;

  fn add(self, c: f64) -> Self::Output {
    Self::Output::new(self.x + c, self.y + c, self.z + c)
  }
}

impl AddAssign for Vector {
  fn add_assign(&mut self, other: Self) {
    *self = *self + &other;
  }
}

impl<'a, 'b> Sub<&'b Vector> for &'a Vector {
  type Output = Vector;

  fn sub(self, other: &'b Vector) -> Self::Output {
    Self::Output::new(self.x - other.x, self.y - other.y, self.z - other.z)
  }
}

impl Mul<f64> for Vector {
  type Output = Self;

  fn mul(self, c: f64) -> Self::Output {
    Self::Output::new(self.x * c, self.y * c, self.z * c)
  }
}

impl Mul<Vector> for f64 {
  type Output = Vector;

  fn mul(self, other: Vector) -> Self::Output {
    other * self
  }
}

impl Div<f64> for Vector {
  type Output = Self;

  fn div(self, c: f64) -> Self::Output {
    Self::Output::new(self.x / c, self.y / c, self.z / c)
  }
}

impl<'a> Div<f64> for &'a Vector {
  type Output = Vector;

  fn div(self, c: f64) -> Self::Output {
    Self::Output::new(self.x / c, self.y / c, self.z / c)
  }
}

impl DivAssign<f64> for Vector {
  fn div_assign(&mut self, c: f64) {
    *self = *self / c
  }
}

impl Sum for Vector {
  fn sum<I>(iter: I) -> Self where I: Iterator<Item = Self> {
    let mut zero = I::Item::zero();

    for vector in iter {
      zero += vector;
    }

    zero
  }
}

impl Neg for Vector {
  type Output = Vector;

  fn neg(self) -> Self::Output {
    Self::Output::new(-self.x, -self.y, -self.z)
  }
}
