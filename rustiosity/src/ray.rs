use vector::Vector;

pub struct Ray {
  pub org: Vector,
  pub dir: Vector,
}

impl Ray {
  pub fn new(org: &Vector, dir: &Vector) -> Ray {
    Ray { org: org.to_owned(), dir: dir.to_owned() }
  }
}
