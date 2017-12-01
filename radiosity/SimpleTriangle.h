#ifndef __SIMPLE_TRIANGLE_H__
#define __SIMPLE_TRIANGLE_H__

#include "Vector.h"
#include "Ray.h"

struct SimpleTriangle {
  Vector a, b, c;
  Vector normal;

  double area;

  SimpleTriangle(const Vector& a_, const Vector& b_, const Vector& c_);
  SimpleTriangle(const Vector &a_, const Vector &b_,
    const Vector &c_, const Color &emission, const Color &color);

  double intersect(const Ray &ray) const;
  const Vector barycentricCoordinatesAt(const Vector &p) const;
};

#endif //  __SIMPLE_TRIANGLE_H__
