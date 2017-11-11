#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <memory>

#include "Shape.h"

struct Triangle: public Shape {
  Vector a, b, c;
  double ab, bc, ca;
  int division_number;

  Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color);

  const unique_ptr<Vector> intersect(const Ray &ray) const;
  void init_patches(const int division_number);
  Color sample_patch(int ia, int ib) const;
};

#endif // __TRIANGLE_H__