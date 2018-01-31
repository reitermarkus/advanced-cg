#include "Sampler.h"

#include <cmath>
#include "macro.h"

// Generates a random sample point on a unit sphere, or a random spherical direction.
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

// Generates a random ray originating from the surface of a sphere.
Ray Sampler::sphericalRay(const Vector &position, double radius) {
  auto sample = pointOnSphere();

  auto point = position + sample * radius;
  auto direction = sample;

  return Ray(point, direction);
}

// Generates a random direction with a given maximum angle in respect to initial direction.
Vector Sampler::randomDirection(Vector direction, double max_angle) {
  // Set up local orthogonal coordinate system u, v, w.
  Vector w = direction;
  Vector u = fabs(w.x) > 0.1 ? Vector(0.0, 1.0, 0.0) : Vector(1.0, 0.0, 0.0);
  u = (u.crossProduct(w)).normalize();
  Vector v = w.crossProduct(u);

  double eps1 = drand48();
  double eps2 = drand48();
  double cos_a = 1.0 - eps1 + eps1 * cos(max_angle);
  double sin_a = sqrt(1.0 - cos_a * cos_a);
  double phi = 2.0 * M_PI * eps2;

  Vector l = u * cos(phi) * sin_a +
             v * sin(phi) * sin_a +
             w * cos_a;

  return l.normalize();
}
