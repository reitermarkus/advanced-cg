#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <memory>

#include "Shape.h"

struct Triangle: public Shape {
  Vector a, b, c;
  Vector b_rel, c_rel;
  double ab, bc, ca;
  int divisions;
  vector<vector<Vector>> subTriangles;

  Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color);

  double intersect(const Ray &ray) const;
  void init_patches(const int divisions);
  vector<vector<Vector>> divide(const int divisions);
  static Vector sample(Vector &p0, Vector &p1, Vector &p2);
};

#endif // __TRIANGLE_H__
