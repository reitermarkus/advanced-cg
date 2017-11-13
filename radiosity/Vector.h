#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <cmath>
#include <iostream>

using namespace std;

typedef struct Vector Color;

struct Vector {
    Vector(const Vector &b);
    Vector(double x_ = 0, double y_ = 0, double z_ = 0);

    Vector operator+(const Vector &b) const;
    Vector operator-(const Vector &b) const;
    Vector operator/(double c) const;
    Vector operator*(double c) const;

    friend Vector operator*(double c, const Vector &b) { return b * c; }

    Vector entrywiseProduct(const Vector &b) const;

    const double lengthSquared() const;
    const double length() const;
    const double distance(const Vector &b) const;
    const Vector normalize() const;
    const double dotProduct(const Vector &b) const;
    const Vector crossProduct(const Vector &b) const;

    double x, y, z;
};

ostream& operator<<(std::ostream &strm, const Vector &v);

#endif //__VECTOR_H__
