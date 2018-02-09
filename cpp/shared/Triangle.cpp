#include "Triangle.h"


static Vector calculateFaceNormal(const Vector a, const Vector b, const Vector c) {
  Vector u = b - a;
  Vector v = c - a;

  return u.crossProduct(v).normalize();
}

Triangle::Triangle(const Vector &a_, const Vector &b_, const Vector &c_, const Color &emission_, const Color &color_, const Refl_t refl_)
          : SceneObject(color_, emission_, refl_), a(a_), b(a_ + b_), c(a_ + c_) {
  this->normal = calculateFaceNormal(a, b, c);
}

double Triangle::intersect(const Ray &ray) const {
  static const double EPSILON = 0.0000001;

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

ostream& operator<<(std::ostream &strm, const Triangle &t) {
  return strm << "Triangle:" << std::endl << t.a << std::endl << t.b <<
    std::endl << t.c << std::endl << t.emission << std::endl << t.color << std::endl;
}
