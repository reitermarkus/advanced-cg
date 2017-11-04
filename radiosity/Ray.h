#ifndef __RAY_H__
#define __RAY_H__

struct Ray {
  Vector org, dir;
  Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_) {}
};

#endif //__RAY_H__
