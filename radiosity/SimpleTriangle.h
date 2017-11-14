#ifndef __SIMPLE_TRIANGLE_H__
#define __SIMPLE_TRIANGLE_H__

#include "Shape.h"

struct SimpleTriangle : public Shape {
  Vector a, b, c;

  SimpleTriangle(const Vector& a_, const Vector& b_, const Vector& c_);
  SimpleTriangle(const Vector &a_, const Vector &b_,
    const Vector &c_, const Color &emission, const Color &color);
};

#endif //  __SIMPLE_TRIANGLE_H__
