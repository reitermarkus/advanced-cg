#ifndef __RAY_H__
#define __RAY_H__

#include "Vector.h"

struct Ray {
  Vector org, dir;
  Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_) {}
};

ostream& operator<<(std::ostream &strm, const Ray &r);

#endif //__RAY_H__
