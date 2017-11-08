#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <cmath>
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

    Vector MultComponents(const Vector &b) const;

    const double LengthSquared() const;
    const double Length() const;
    const double distance(const Vector &b) const;
    const Vector Normalized() const;
    const double Dot(const Vector &b) const;
    const Vector Cross(const Vector &b) const;

    double x, y, z;
};

#endif //__VECTOR_H__
