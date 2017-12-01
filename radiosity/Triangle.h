#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <memory>
#include <vector>

#include "PatchTriangle.h"
#include "SimpleTriangle.h"

struct Triangle: public SimpleTriangle {
  Vector b_rel, c_rel;

  Color emission, color;

  int divisions;
  vector<PatchTriangle> subTriangles;
  vector<Color> patch;

  Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color);

  void init_patches(const int divisions);
  vector<PatchTriangle> divide(const int divisions);
  static Vector sample(Vector &p0, Vector &p1, Vector &p2);
};

#endif // __TRIANGLE_H__
