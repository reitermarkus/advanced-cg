#include "SimpleTriangle.h"

#include <limits>

static double calculateArea(const Vector a, const Vector b, const Vector c) {
  return ((b - a).crossProduct(c - a) / 2.0).length();
}

static Vector calculateFaceNormal(const Vector a, const Vector b, const Vector c) {
  Vector u = b - a;
  Vector v = c - a;

  return u.crossProduct(v);
}

SimpleTriangle::SimpleTriangle(const Vector& a_, const Vector& b_, const Vector& c_) : a(a_), b(b_), c(c_) {
  normal = calculateFaceNormal(a, b, c).normalize();
  area = calculateArea(a, b, c);
}

//  Tests whether a Ray intersects with a Triangle using the Möller-Trumbore algorithm.
//  Either returns 0.0 or the distance on the Ray.
//
double SimpleTriangle::intersect(const Ray &ray) const {
  static const double EPSILON = numeric_limits<double>::epsilon();

  Vector edge_1 = this->b - this->a;
  Vector edge_2 = this->c - this->a;

  Vector h = ray.dir.crossProduct(edge_2);
  double a = edge_1.dotProduct(h);

  if (a > -EPSILON && a < EPSILON)
    return 0.0;

  double f = 1.0 / a;
  Vector s = ray.org - this->a;
  double u = f * s.dotProduct(h);

  if (u < 0.0 || u > 1.0)
    return 0.0;

  Vector q = s.crossProduct(edge_1);
  double v = f * ray.dir.dotProduct(q);

  if (v < 0.0 || u + v > 1.0)
    return 0.0;

  double t = f * edge_2.dotProduct(q);

  if (t <= EPSILON)
    return 0.0;

  return t;
}

const Vector SimpleTriangle::barycentricCoordinatesAt(const Vector &p) const {
  double area_p_b_c = ((b - p).crossProduct(c - p) / 2.0).length();
  double area_p_b_a = ((b - a).crossProduct(p - a) / 2.0).length();

  double alpha = area_p_b_c / area;
  double beta  = area_p_b_a / area;
  double gamma = 1.0 - alpha - beta;

  return Vector(alpha, beta, gamma);
}
