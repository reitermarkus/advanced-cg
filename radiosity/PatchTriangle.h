#ifndef __PATCH_TRIANGLE_H__
#define __PATCH_TRIANGLE_H__

#include "Vector.h"
#include "SimpleTriangle.h"

struct PatchTriangle : public SimpleTriangle {
  Color color;
  PatchTriangle(const Vector& a_, const Vector& b_, const Vector& c_) : SimpleTriangle(a_, b_, c_) {}
};

#endif // __PATCH_TRIANGLE_H__
