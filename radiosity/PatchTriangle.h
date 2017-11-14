#ifndef __PATCH_TRIANGLE_H__
#define  __PATCH_TRIANGLE_H__

#include "Vector.h"

struct PatchTriangle {
  Vector a, b, c;
  PatchTriangle(const Vector& a_, const Vector& b_, const Vector& c_) : a(a_), b(b_), c(c_) {}
};

#endif // __PATCH_TRIANGLE_H__