#ifndef __SHAPE_H__
#define __SHAPE_H__

#include <vector>

#include "Vector.h"
#include "Ray.h"

struct Shape {
  Color emission, color;
  Vector normal;
  double area;

  vector<Color> patch;

  Shape(const Color &emission_, const Color &color_);
  Shape() = default;
};

#endif // __SHAPE_H__
