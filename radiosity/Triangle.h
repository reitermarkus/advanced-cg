#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <memory>

#include "Shape.h"

struct Triangle: public Shape {
  Vector a, b, c;
  Vector b_rel, c_rel;
  double ab, bc, ca;
  int division_number;

  Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color);

  double intersect(const Ray &ray) const;
  void init_patches(const int division_number);
  unique_ptr<vector<vector<Vector>>> divide(int divisions) const;
};

#endif // __TRIANGLE_H__
