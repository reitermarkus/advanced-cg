#include "Triangle.h"

#include <cmath>

#include "Ray.h"
#include "Vector.h"

static Vector calculateFaceNormal(const Vector a, const Vector b, const Vector c) {
  Vector u = b - a;
  Vector v = c - a;

  return u.crossProduct(v);
}

Triangle::Triangle(const Vector &a_, const Vector &b_, const Vector &c_,
           const Color &emission, const Color &color): Shape(emission, color) {
  a = a_;
  b = a_ + b_; // b_ is given relative to a_
  c = a_ + c_; // c_ is given relative to a_
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

//  Tests whether a Ray intersects with a Triangle using the Möller-Trumbore algorithm.
//  Either returns a pointer to the intersection Vector or null.
//
const unique_ptr<Vector> Triangle::intersect(const Ray &ray) const {
  static const double EPSILON = 0.0000001;

  Vector edge_1 = this->b - this->a;
  Vector edge_2 = this->c - this->a;

  Vector h = ray.dir.crossProduct(edge_2);
  double a = edge_1.dotProduct(h);

  if (a > -EPSILON && a < EPSILON)
    return nullptr;

  double f = 1 / a;
  Vector s = ray.org;
  double u = f * s.dotProduct(h);

  if (u < 0.0 || u > 1.0)
    return nullptr;

  Vector q = s.crossProduct(edge_1);
  double v = f * ray.dir.dotProduct(q);

  if (v < 0.0 || u + v > 1.0)
    return nullptr;

  double t = f * edge_2.dotProduct(q);

  if (t <= EPSILON)
    return nullptr;

  return make_unique<Vector>(ray.org + ray.dir * t);
}

void Triangle::init_patches(const int division_number) {
  this->division_number = division_number;
  this->patch.clear();
  this->patch.resize(pow(division_number, 2));
}

void Triangle::divide(int divisions) const {
  auto delta_x = b_rel / divisions;
  auto delta_y = c_rel / divisions;

  unique_ptr<vector<vector<Vector>>> subTriangles = make_unique<vector<vector<Vector>>>();

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

      subTriangles->push_back({v1, v2, v3});

      cout << v1 << ", " << v2 << ", " << v3 << ", " << endl;
    }

    row++;

    if (row % 2 == 1) {
      row_size--;
    }
  }

  cout << subTriangles->size() << endl;
}

Color Triangle::sample_patch(int ia, int ib) const {
  return this->patch[0];
}
