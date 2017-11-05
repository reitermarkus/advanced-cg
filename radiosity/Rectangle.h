#ifndef __RECTANGLE_H__
#define __RECTANGLE_H__

#include "Shape.h"

struct Rectangle: public Shape {
  Rectangle(const Vector p0_, const Vector &a_, const Vector &b_,
            const Color &emission_, const Color &color_):
    Shape(p0_, a_, b_, emission_, color_) {}

  const double intersect(const Ray &ray) const;
};

#endif // __RECTANGLE_H__
