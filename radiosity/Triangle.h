#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <memory>

#include "Shape.h"

struct Triangle: public Shape {
  Vector a, b, c;

  Triangle(const Vector a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color):
    Shape(emission, color), a(a_), b(b_), c(c_) {}

  const unique_ptr<Vector> intersect(const Ray &ray) const;
};

#endif // __TRIANGLE_H__
