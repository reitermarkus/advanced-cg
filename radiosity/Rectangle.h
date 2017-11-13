#ifndef __RECTANGLE_H__
#define __RECTANGLE_H__

#include "Shape.h"

struct Rectangle: public Shape {
  Vector p0;
  Vector edge_a, edge_b;

  int a_num, b_num;
  double a_len, b_len;

  Rectangle(const Vector p0_, const Vector &a_, const Vector &b_,
            const Color &emission_, const Color &color_);

  void init_patches(const int a_num_, const int b_num_);
  Color sample_patch(int ia, int ib) const;
  double intersect(const Ray &ray) const;
};

#endif // __RECTANGLE_H__
