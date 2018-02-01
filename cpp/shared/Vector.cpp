#include "Vector.h"

#include <algorithm>

Vector::Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
Vector::Vector(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

Vector& Vector::operator+=(const Vector &b) {
  x += b.x;
  y += b.y;
  z += b.z;

  return *this;
}

Vector Vector::operator-(const Vector &b) const {
  return Vector(x - b.x, y - b.y, z - b.z);
}

Vector Vector::operator-() const {
  return Vector(-x, -y, -z);
}

Vector& Vector::operator/=(double c) {
  x /= c;
  y /= c;
  z /= c;

  return *this;
}

Vector Vector::operator*(double c) const {
  return Vector(x * c, y * c, z * c);
}

bool Vector::operator==(const Vector &other) const {
  return x == other.x && y == other.y && z == other.z;
}

bool Vector::operator<(const Vector &other) const {
  if (x < other.x) {
    return true;
  } else if (y < other.y) {
    return true;
  } else if (z < other.z) {
    return true;
  }

  return false;
}

Vector Vector::entrywiseProduct(const Vector &b) const {
  return Vector(x * b.x, y * b.y, z * b.z);
}

double Vector::lengthSquared() const {
  return dotProduct(*this);
}

double Vector::length() const {
  return sqrt(lengthSquared());
}

double Vector::distance(const Vector &b) const {
  return (*this - b).length();
}

Vector Vector::normalize() const {
  return Vector(x, y, z) / length();
}

double Vector::dotProduct(const Vector &b) const {
  return x * b.x + y * b.y + z * b.z;
}

Vector Vector::crossProduct(const Vector &b) const {
  return Vector(
    (y * b.z) - (z * b.y),
    (z * b.x) - (x * b.z),
    (x * b.y) - (y * b.x)
  );
}

Vector Vector::clamp(double min, double max) {
  return Vector(
    std::clamp(x, min, max),
    std::clamp(y, min, max),
    std::clamp(z, min, max)
  );
}

double Vector::max() const {
  return std::max(x, std::max(y, z));
}

double Vector::min() const {
  return std::min(x, std::min(y, z));
}

ostream& operator<<(std::ostream &strm, const Vector &v) {
  return strm << "Vector(" << v.x << ", " << v.y << ", " << v.z << ")";
}
