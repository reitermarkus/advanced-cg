#ifndef __SCENE_OBJECT_H__
#define __SCENE_OBJECT_H__

#include "../shared/Vector.h"
#include "../shared/Ray.h"

struct SceneObject {
  Color color, emission;
  const bool isSphere;

  SceneObject(const Vector color_, const Vector emission_, const bool isSphere_ = false) : color(color_), emission(emission_), isSphere(isSphere_) {}

  virtual double intersect(const Ray &ray) const = 0;

};

#endif // __SCENE_OBJECT_H__
