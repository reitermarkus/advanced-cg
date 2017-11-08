#include "Vector.h"

Vector::Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
Vector::Vector(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

Vector Vector::operator+(const Vector &b) const {
  return Vector(x + b.x, y + b.y, z + b.z);
}

Vector Vector::operator-(const Vector &b) const {
  return Vector(x - b.x, y - b.y, z - b.z);
}

Vector Vector::operator/(double c) const {
  return Vector(x / c, y / c, z / c);
}

Vector Vector::operator*(double c) const {
  return Vector(x * c, y * c, z * c);
}

Vector Vector::MultComponents(const Vector &b) const {
  return Vector(x * b.x, y * b.y, z * b.z);
}

const double Vector::LengthSquared() const {
  return x * x + y * y + z * z;
}

const double Vector::Length() const {
  return sqrt(LengthSquared());
}

const double Vector::distance(const Vector &b) const {
  return sqrt(pow(x - b.x, 2) + pow(y - b.y, 2) + pow(z - b.z, 2));
}

const Vector Vector::Normalized() const {
  return Vector(x, y, z) / sqrt(x * x + y * y + z * z);
}

const double Vector::Dot(const Vector &b) const {
  return x * b.x + y * b.y + z * b.z;
}

const Vector Vector::Cross(const Vector &b) const {
  return Vector((y * b.z) - (z * b.y), (z * b.x) - (x * b.z),
                (x * b.y) - (y * b.x));
}
