#ifndef __SHAPE_H__
#define __SHAPE_H__

#include <vector>

#include "Vector.h"
#include "Ray.h"

struct Shape {
  Color emission, color;
  Vector normal;

  vector<Color> patch;

  Shape(const Color &emission_, const Color &color_);

  virtual void init_patchs(const int a_num_, const int b_num_) = 0;
  virtual Color sample_patch(int ia, int ib) const = 0;
};

#endif // __SHAPE_H__
