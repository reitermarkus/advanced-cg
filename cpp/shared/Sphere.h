#ifndef __SPHERE_H__
#define __SPHERE_H__

#include "../shared/Vector.h"
#include "../shared/Ray.h"

#include "SceneObject.h"

struct Sphere : public SceneObject {
  double radius;
  Vector position;

  Sphere(double radius_, const  Vector &position_, const Vector &emission_, const Vector &color_, const Refl_t refl_) :
    SceneObject(color_, emission_, refl_, true),
    radius(radius_), position(position_) {}

  double intersect(const Ray &ray) const override;
};

#endif // __SPHERE_H__
