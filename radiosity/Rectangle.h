#ifndef __RECTANGLE_H__
#define __RECTANGLE_H__

#include <vector>

#include "Vector.h"
#include "Ray.h"

struct Rectangle {
  Vector p0;
  Vector edge_a, edge_b;
  Color emission, color;
  Vector normal;

  vector<Color> patch;
  int a_num, b_num;
  double a_len, b_len;

  Rectangle(const Vector p0_, const Vector &a_, const Vector &b_,
            const Color &emission_, const Color &color_);

  Color sample_patch(int ia, int ib) const;
  void init_patchs(const int a_num_, const int b_num_);
  const double intersect(const Ray &ray);
};

#endif //__RECTANGLE_H__