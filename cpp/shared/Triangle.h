#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "SceneObject.h"

struct Triangle: public SceneObject {
  Vector a, b, c;
  Vector normal;

  Triangle(const Vector &a_, const Vector &b_, const Vector &c_, const Color &emission_, const Color &color_, const Refl_t refl_);

  double intersect(const Ray &ray) const override;
};

ostream& operator<<(std::ostream &strm, const Triangle &t);

#endif // __TRIANGLE_H__
