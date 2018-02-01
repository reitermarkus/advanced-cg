#include <iostream>
#include <vector>

#include "Sampler.h"
#include "HSV.h"
#include "SceneObject.h"
#include "Triangle.h"
#include "Sphere.h"
#include "Vector.h"
#include "Diffraction.h"
#include "Image.h"

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

Vector perfectReflection(Vector dir, Vector normal) {
  return dir - normal * 2 * normal.dotProduct(dir);
}

pair<long, long> traceRay(Ray &ray, HSV &emission, double aperture, double focal_length) {
  double t;
  size_t id = 0;

  if (!intersect(ray, t, id)) { // No intersection with scene.
    return make_pair(-1, -1);
  }

  const SceneObject* obj = objects[id];

  Vector hitpoint = ray.org + ray.dir * t;    /* Intersection point */

  /* Normal at intersection */
  Vector normal;

  if (obj->isSphere) {
    const Sphere* sphere = static_cast<const Sphere*>(obj);
    normal = (hitpoint - sphere->position).normalize();
  } else {
    const Triangle* tri = static_cast<const Triangle*>(obj);
    normal = tri->normal;
  }

  Vector nl = normal;

  /* Obtain flipped normal, if object hit from inside */
  if (normal.dotProduct(ray.dir) >= 0)
    nl = -nl;

  Ray reflection_ray = Ray(hitpoint, perfectReflection(ray.dir, normal));

  auto wavelength = emission.hueAsWavelength();
  auto min_diffraction_angle = minDiffractionAngle(wavelength, aperture, focal_length);

  return make_pair(0, 0);
}


int main() {
  srand(time(0));

  auto width = 1024;
  auto height = 768;

  double aperture = 2.6;
  double focal_length = 120.0;

  auto samples = 4;

  for (auto &light : lights) {
    objects.push_back(&light);
  }

  vector<Triangle> layers;


  for(int i = 0; i < 5; i++) {
    const auto h = Vector( 0.0, 100.0,  0.0);
    const auto w = Vector(100.0,  0.0,  0.0);
    const auto d = Vector( 0.0,  0.0, 1.0);

    const auto color = Color(1.0, 1.0, 1.0);

    const auto position = Vector(0.0, 0.0, 100 + i);

    layers.push_back(Triangle(i * d + position, w + h,     h, Color(), color, REFR)); // Front
    layers.push_back(Triangle(i * d + position,     w, w + h, Color(), color, REFR)); // Front

    layers.push_back(Triangle(i * d + d + position + w,     -w, -w + h, Color(), color, REFR)); // Back
    layers.push_back(Triangle(i * d + d + position + w, -w + h,      h, Color(), color, REFR)); // Back
  }

  for(const auto& layer : layers) {
    objects.push_back(&layer);
  }


  Image img(width, height);

  long hits = 0;

  auto ray_samples = width * height * samples;

  for (auto &light : lights) {
    auto emission = HSV::rgbToHsv(Color(light.emission.x, light.emission.y, light.emission.z));

    #pragma omp parallel for schedule(dynamic, 1)
    for (auto s = 0; s < ray_samples; s++) {
      auto ray = Sampler::sphericalRay(light.position, light.radius);
      auto ray_emission = HSV::withRandomHue(emission.s, emission.v);

      auto pixels = traceRay(ray, ray_emission, aperture, focal_length);

      if (!(pixels.first == -1 || pixels.second == -1)) {
        #pragma omp critical
        {
          cout << pixels.first << ", " << pixels.second << endl;
          img.addColor(pixels.first, pixels.second, light.emission);
          hits++;
        }


      }
    }
  }

  cout << hits << "/" << ray_samples << endl;
  cout << (double)hits / (double)ray_samples * 100 << " %" << endl;

  img.save("image");
}
