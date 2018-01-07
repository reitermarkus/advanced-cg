#ifndef __TRIGONOMETRICAL_H__
#define __TRIGONOMETRICAL_H__

#include "../shared/Vector.h"
#include "../shared/Ray.h"

struct Trigonometrical {
  Color color, emission;

  Trigonometrical(const Vector color_, const Vector emission_) : color(color_), emission(emission_) { }

  virtual double intersect(const Ray &ray) = 0;
};

#endif // __TRIGONOMATRICAL_H__
