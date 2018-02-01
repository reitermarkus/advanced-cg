#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <cmath>
#include <iostream>

using namespace std;

typedef struct Vector Color;

struct Vector {
    Vector(const Vector &b);
    Vector(double x_ = 0, double y_ = 0, double z_ = 0);

    Vector& operator+=(const Vector &b);
    Vector operator-(const Vector &b) const;
    Vector operator-() const;
    Vector& operator/=(double c);
    Vector operator*(double c) const;
    bool operator==(const Vector &b) const;
    bool operator<(const Vector &b) const;

    friend Vector operator+(Vector a, const Vector &b) { return a += b; }
    friend Vector operator*(double c, const Vector &b) { return b * c; }
    friend Vector operator/(Vector a, double c) { return a /= c; }

    double lengthSquared() const;
    double length() const;
    double distance(const Vector &b) const;
    Vector normalize() const;
    double dotProduct(const Vector &b) const;
    Vector crossProduct(const Vector &b) const;
    Vector entrywiseProduct(const Vector &b) const;
    Vector clamp(double min, double max);
    double max() const;
    double min() const;

    double x, y, z;
};

ostream& operator<<(std::ostream &strm, const Vector &v);

namespace std {
  template<>
  struct hash<Vector> {
    size_t operator()(const Vector& v) const {
      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:
      return ((hash<double>()(v.x) ^ (hash<double>()(v.y) << 1)) >> 1) ^ (hash<double>()(v.z) << 1);
    }
  };
}

#endif //__VECTOR_H__
