#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "Shape.h"

struct Triangle: public Shape {
  Triangle(const Vector p0, const Vector &edge_a, const Vector &edge_b,
           const Color &emission, const Color &color):
    Shape(p0, edge_a, edge_b, emission, color) {}

  const double intersect(const Ray &ray) const;
};

#endif // __TRIANGLE_H__
