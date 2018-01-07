#include "Triangle.h"

Triangle::Triangle(const Vector &a_, const Vector &b_,
          const Vector &c_, const Color &emission_, const Color &color_)
          : SceneObject(color_, emission_), a(a_), b(b_), c(c_) {}

double Triangle::intersect(const Ray &ray) const {
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
