#include "Triangle.h"

#include <cmath>

#include "Ray.h"
#include "Vector.h"

static Vector calculateFaceNormal(const Vector a, const Vector b, const Vector c) {
  Vector u = b - a;
  Vector v = c - a;

  return u.crossProduct(v);
}

Triangle::Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color): Shape(emission, color), a(a_), b(b_), c(c_) {
  normal = calculateFaceNormal(a_, b_, c_).normalize();

  ab = a.distance(b);
  bc = b.distance(c);
  ca = c.distance(a);

  // Calculate area of triangle using Heron's Formula.
  auto s = (ab + bc + ca) / 2.0;
  area = sqrt(s * (s - ab) * (s - bc) * (s - ca));
}

//  Tests whether a Ray intersects with a Triangle using the Möller-Trumbore algorithm.
//  Either returns a pointer to the intersection Vector or null.
//
const unique_ptr<Vector> Triangle::intersect(const Ray &ray) const {
  static const double EPSILON = 0.0000001;

  Vector edge_1 = this->b - this->a;
  Vector edge_2 = this->c - this->a;

  Vector h = ray.dir.crossProduct(edge_2);
  double a = edge_1.dotProduct(h);

  if (a > -EPSILON && a < EPSILON)
    return nullptr;

  double f = 1 / a;
  Vector s = ray.org;
  double u = f * s.dotProduct(h);

  if (u < 0.0 || u > 1.0)
    return nullptr;

  Vector q = s.crossProduct(edge_1);
  double v = f * ray.dir.dotProduct(q);

  if (v < 0.0 || u + v > 1.0)
    return nullptr;

  double t = f * edge_2.dotProduct(q);

  if (t <= EPSILON)
    return nullptr;

  return make_unique<Vector>(ray.org + ray.dir * t);
}
