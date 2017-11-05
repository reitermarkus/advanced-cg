#include "Shape.h"

Shape::Shape(const Vector p0_, const Vector &a_, const Vector &b_, const Color &emission_, const Color &color_) {
  p0 = p0_;
  edge_a = a_;
  edge_b = b_;
  emission = emission_;
  color = color_;
  normal = edge_a.Cross(edge_b);
  normal = normal.Normalized();
  a_len = edge_a.Length();
  b_len = edge_b.Length();
}

void Shape::init_patchs(const int a_num_, const int b_num_) {
  a_num = a_num_;
  b_num = b_num_;
  patch.clear();
  patch.resize(a_num * b_num);
}

Color Shape::sample_patch(int ia, int ib) const {
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
