#ifndef __SPHERE_H__
#define __SPHERE_H__

#include "../shared/Vector.h"
#include "../shared/Ray.h"

/*------------------------------------------------------------------
| Scene objects are spheres; material either perfectly diffuse,
| specular (mirror reflection) or transparent (refraction/reflection)
| (DIFFuse, SPECular, REFRactive)
------------------------------------------------------------------*/
enum Refl_t { DIFF, SPEC, REFR };

struct Sphere {
    double radius;
    Vector position;
    Color emission, color;
    Refl_t refl;

    Sphere(double radius_, Vector position_, Vector emission_,
           Vector color_, Refl_t refl_):
           radius(radius_), position(position_), emission(emission_),
           color(color_), refl(refl_) {}

    double Intersect(const Ray &ray) const;
};

#endif // __SPHERE_H__
