#ifndef __SHAPE_H__
#define __SHAPE_H__

#include <vector>

#include "Vector.h"
#include "Ray.h"

struct Shape {
  Vector p0;
  Vector edge_a, edge_b;
  Color emission, color;
  Vector normal;

  vector<Color> patch;
  int a_num, b_num;
  double a_len, b_len;

  Shape(const Vector p0_, const Vector &a_, const Vector &b_,
        const Color &emission_, const Color &color_);

  void init_patchs(const int a_num_, const int b_num_);
  Color sample_patch(int ia, int ib) const;
  virtual const double intersect(const Ray &ray) const = 0;
};

#endif // __SHAPE_H__
