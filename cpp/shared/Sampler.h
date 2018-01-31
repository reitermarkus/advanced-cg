#ifndef __SAMPLER_H__
#define __SAMPLER_H__

#include "Vector.h"

struct Sampler {
  static Vector pointOnSphere();
  static Vector randomDirection(Vector direction, double max_angle);
};

#endif // __SAMPLER_H__
