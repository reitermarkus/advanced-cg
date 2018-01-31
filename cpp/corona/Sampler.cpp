#include "Sampler.h"

#include <cmath>
#include "macro.h"

// http://mathworld.wolfram.com/SpherePointPicking.html
Vector Sampler::pointOnSphere() {
  auto u = drand48();
  auto v = drand48();

  auto theta = 2 * M_PI * u;
  auto phi = acos(2 * v - 1.0);

  u = cos(phi);

  auto sqrt_1_u_2 = sqrt(1.0 - pow(u, 2));

  auto x = sqrt_1_u_2 * cos(theta);
  auto y = sqrt_1_u_2 * sin(theta);
  auto z = u;

  return Vector(x, y, z);
}
