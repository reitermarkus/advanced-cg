#ifndef __SAMPLER_H__
#define __SAMPLER_H__

#include "Vector.h"
#include "Ray.h"

struct Sampler {
  static Vector pointOnSphere();
  static Ray sphericalRay(const Vector &position, double radius);
  static Vector randomDirection(Vector direction, double max_angle);
};

#endif // __SAMPLER_H__
