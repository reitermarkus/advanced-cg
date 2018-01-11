#ifndef __SCENE_OBJECT_H__
#define __SCENE_OBJECT_H__

#include "../shared/Vector.h"
#include "../shared/Ray.h"

/*------------------------------------------------------------------
| Scene objects are spheres or triangles; material either perfectly diffuse,
| specular (mirror reflection) or transparent (refraction/reflection)
| (DIFFuse, SPECular, REFRactive)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR, GLOS, TRAN };

struct SceneObject {
  Color color, emission;
  Refl_t refl;

  const bool isSphere;

  SceneObject(const Vector &color_, const Vector &emission_, const Refl_t refl_, const bool isSphere_ = false) : color(color_), emission(emission_), refl(refl_), isSphere(isSphere_) {}

  virtual double intersect(const Ray &ray) const = 0;

};

#endif // __SCENE_OBJECT_H__
