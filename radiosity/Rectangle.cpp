#include "Rectangle.h"

const double Rectangle::intersect(const Ray &ray) const {
  const double t = (p0 - ray.org).Dot(normal) / ray.dir.Dot(normal);
  if (t <= 0.00001)
    return 0.0;

  Vector p = ray.org + ray.dir * t;
  Vector d = p - p0;
  const double ddota = d.Dot(edge_a);
  if (ddota < 0.0 || ddota > edge_a.LengthSquared())
    return 0.0;

  const double ddotb = d.Dot(edge_b);
  if (ddotb < 0.0 || ddotb > edge_b.LengthSquared())
    return 0.0;

  return t;
}
