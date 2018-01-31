#include <iostream>
#include <vector>

#include "Sampler.h"
#include "HSV.h"
#include "SceneObject.h"
#include "Sphere.h"
#include "Vector.h"

using namespace std;

vector<const SceneObject*> objects = vector<const SceneObject*>();

vector<Sphere> lights = {
  Sphere(4.0, Vector(50, 81.6 - 16.5, 81.6), Color(4, 4, 4) * 100, Color(), DIFF), /* Light */
};


/******************************************************************
* Check for closest intersection of a ray with the scene;
* returns true if intersection is found, as well as ray parameter
* of intersection and id of intersected object
*******************************************************************/
bool intersect(const Ray &ray, double &t, size_t &id) {
  t = 1e20;

  for (size_t i = 0; i < objects.size(); i++) {
    double d = objects[i]->intersect(ray);
    if (d > 0.0  && d < t) {
      t = d;
      id = i;
    }
  }

  return t < 1e20;
}

void traceRay(Ray &ray, HSV &emission) {
  double t;
  size_t id = 0;

  if (!intersect(ray, t, id)) { // No intersection with scene.
    // return Color(0.0, 0.0, 0.0);
  }

  auto wavelength = emission.hueAsWavelength();
}


int main() {
  srand(time(0));

  auto width = 1024;
  auto height = 768;

  auto samples = 4;

  for (auto &light : lights) {
    objects.push_back(&light);
  }

  auto ray_samples = width * height * samples;

  for (auto &light : lights) {
    auto emission = HSV::from(RGB(light.emission.x, light.emission.y, light.emission.z));

    for (auto s = 0; s < ray_samples; s++) {
      auto ray = Sampler::sphericalRay(light.position, light.radius);
      auto ray_emission = HSV::withRandomHue(emission.s, emission.v);

      traceRay(ray, ray_emission);
    }
  }

  return 0;
}
