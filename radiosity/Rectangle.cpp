#include "Rectangle.h"

Rectangle::Rectangle(const Vector p0_, const Vector &a_, const Vector &b_,
          const Color &emission_, const Color &color_)
    : p0(p0_), edge_a(a_), edge_b(b_), emission(emission_), color(color_) {
  normal = edge_a.Cross(edge_b);
  normal = normal.Normalized();
  a_len = edge_a.Length();
  b_len = edge_b.Length();
}

Color Rectangle::sample_patch(int ia, int ib) const {
  if (ia < 0)
    ia = 0;
  if (ia >= a_num)
    ia = a_num - 1;
  if (ib < 0)
    ib = 0;
  if (ib >= b_num)
    ib = b_num - 1;
  return patch[ia * b_num + ib];
}

void Rectangle::init_patchs(const int a_num_, const int b_num_) {
  a_num = a_num_;
  b_num = b_num_;
  patch.clear();
  patch.resize(a_num * b_num);
}

const double Rectangle::intersect(const Ray &ray) {
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