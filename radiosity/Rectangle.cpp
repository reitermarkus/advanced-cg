#include "Rectangle.h"

Rectangle::Rectangle(const Vector p0_, const Vector &a_, const Vector &b_,
          const Color &emission_, const Color &color_): Shape(emission_, color_), p0(p0_), edge_a(a_), edge_b(b_) {
  normal = edge_a.crossProduct(edge_b);
  normal = normal.normalize();
  a_len = edge_a.length();
  b_len = edge_b.length();
  area = a_len * b_len;
}

void Rectangle::init_patches(const int a_num_, const int b_num_) {
  a_num = a_num_;
  b_num = b_num_;
  patch.clear();
  patch.resize(a_num * b_num);
}

Color Rectangle::sample_patch(int ia, int ib) const {
  ia = clamp(ia, 0, a_num - 1);
  ib = clamp(ib, 0, b_num - 1);
  return patch[ia * b_num + ib];
}

double Rectangle::intersect(const Ray &ray) const {
  const double t = (p0 - ray.org).dotProduct(normal) / ray.dir.dotProduct(normal);
  if (t <= 0.00001)
    return 0.0;

  Vector p = ray.org + ray.dir * t;
  Vector d = p - p0;
  const double ddota = d.dotProduct(edge_a);
  if (ddota < 0.0 || ddota > edge_a.lengthSquared())
    return 0.0;

  const double ddotb = d.dotProduct(edge_b);
  if (ddotb < 0.0 || ddotb > edge_b.lengthSquared())
    return 0.0;

  return t;
}
