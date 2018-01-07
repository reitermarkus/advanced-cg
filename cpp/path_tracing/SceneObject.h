#ifndef __SCENE_OBJECT_H__
#define __SCENE_OBJECT_H__

#include "../shared/Vector.h"
#include "../shared/Ray.h"

struct SceneObject {
  Color color, emission;

  SceneObject(const Vector color_, const Vector emission_) : color(color_), emission(emission_) { }

  virtual double intersect(const Ray &ray) = 0;
};

#endif // __SCENE_OBJECT_H__
