#ifndef __SPHERE_H__
#define __SPHERE_H__

#include "../shared/Vector.h"
#include "../shared/Ray.h"

#include "Trigonometrical.h"

/*------------------------------------------------------------------
| Scene objects are spheres; material either perfectly diffuse,
| specular (mirror reflection) or transparent (refraction/reflection)
| (DIFFuse, SPECular, REFRactive)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR };

struct Sphere : public Trigonometrical {
  double radius;
  Vector position;
  Refl_t refl;

  Sphere(double radius_, Vector position_, Vector emission_,
          Vector color_, Refl_t refl_);

  double intersect(const Ray &ray) override;
};

#endif // __SPHERE_H__
