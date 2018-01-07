#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "Trigonometrical.h"

struct Triangle: public Trigonometrical {
  Vector a, b, c;

  Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color);

  double intersect(const Ray &ray) override;
};

#endif // __TRIANGLE_H__
