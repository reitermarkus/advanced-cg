#include "Triangle.h"

#include <cmath>

#include "macro.h"

#include "Ray.h"
#include "Vector.h"

static Vector calculateFaceNormal(const Vector a, const Vector b, const Vector c) {
  Vector u = b - a;
  Vector v = c - a;

  return u.crossProduct(v);
}

Triangle::Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color): SimpleTriangle(a_, b_, c_, emission, color) {
  b_rel = b_;
  c_rel = c_;

  normal = calculateFaceNormal(a, b, c).normalize();

  ab = a.distance(b);
  bc = b.distance(c);
  ca = c.distance(a);

  // Calculate area of triangle using Heron's Formula.
  auto s = (ab + bc + ca) / 2.0;
  area = sqrt(s * (s - ab) * (s - bc) * (s - ca));
}

void Triangle::init_patches(const int divisions) {
  this->divisions = divisions;
  this->patch.clear();
  this->patch.resize(pow(divisions, 2));
  this->divide(divisions);
}

vector<PatchTriangle> Triangle::divide(const int divisions) {
  auto delta_x = b_rel / divisions;
  auto delta_y = c_rel / divisions;

  this->subTriangles.clear();

  // Loop through patches, from bottom to top, left to right.
  auto row_size = divisions;
  auto row = 0;
  while (row_size > 0) {
    for (auto x = 0; x < row_size; x++) {
      auto y = int(row / 2);

      auto offset = row % 2;

      auto v1 = a + (x + offset) * delta_x + y * delta_y;
      auto v2 = a + (x + 1) * delta_x + (y + offset) * delta_y;
      auto v3 = a + x * delta_x + (y + 1) * delta_y;

      auto patch_triangle = PatchTriangle(v1, v2, v3);
      this->subTriangles.push_back(patch_triangle);
    }

    row++;

    if (row % 2 == 1) {
      row_size--;
    }
  }

  return this->subTriangles;
}

Vector Triangle::sample(Vector &p0, Vector &p1, Vector &p2) {
  srand(time(nullptr));

  auto epsilon_0 = drand48();
  auto epsilon_1 = drand48();

  auto delta_0 = 1.0 - sqrt(epsilon_0);
  auto delta_1 = epsilon_1 * sqrt(epsilon_0);
  auto delta_2 = 1 - delta_0 - delta_1; // = sqrt(epsilon_0) * (1.0 - epsilon_1)

  Vector q = delta_0 * p0 + delta_1 * p1 + delta_2 * p2;

  return q;
}
