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

Vector Vector::entrywiseProduct(const Vector &b) const {
  return Vector(x * b.x, y * b.y, z * b.z);
}

const double Vector::lengthSquared() const {
  return dotProduct(*this);
}

const double Vector::length() const {
  return sqrt(lengthSquared());
}

const double Vector::distance(const Vector &b) const {
  return (*this - b).length();
}

const Vector Vector::normalize() const {
  return Vector(x, y, z) / length();
}

const double Vector::dotProduct(const Vector &b) const {
  return x * b.x + y * b.y + z * b.z;
}

const Vector Vector::crossProduct(const Vector &b) const {
  return Vector(
    (y * b.z) - (z * b.y),
    (z * b.x) - (x * b.z),
    (x * b.y) - (y * b.x)
  );
}
