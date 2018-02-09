#include <iostream>
#include <vector>
#include <iomanip>

#include "macro.h"

#include "Sampler.h"
#include "HSV.h"
#include "SceneObject.h"
#include "Triangle.h"
#include "Sphere.h"
#include "Vector.h"
#include "Diffraction.h"
#include "Image.h"
#include "refractive_index.h"

using namespace std;

enum Wave { R, G, B };

vector<const SceneObject*> objects = vector<const SceneObject*>();

const auto light = Sphere(10.0, Vector(50.0, 50.0, 0.0), Color(4, 4, 4) * 100, Color(), DIFF);

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

Color radiance(const Ray &ray, Wave wave) {
  double t;
  size_t id = 0;

  if (!intersect(ray, t, id))   /* No intersection with scene */
    return Color(0.0, 0.0, 0.0);

  const SceneObject* obj = objects[id];

  Vector hitpoint = ray.org + ray.dir * t;    /* Intersection point */

  /* Normal at intersection */
  Vector normal;

  if (obj->isSphere) {
    const Sphere* sphere = static_cast<const Sphere*>(obj);
    normal = (hitpoint - sphere->position).normalize();

    // Only the light is a sphere, so return if the ray reaches the light.
    return sphere->emission;
  } else {
    const Triangle* tri = static_cast<const Triangle*>(obj);
    normal = tri->normal;
  }

  Vector nl = normal;

  /* Obtain flipped normal, if object hit from inside */
  if (normal.dotProduct(ray.dir) >= 0)
    nl = -nl;

  double n1 = 1.0;
  double n2 = 1.0;

  auto offset = (drand48() - 0.5) / 10.0;

  switch (wave) {
    case R: n2 = n2 + 0.01 + offset; break;
    case G: n2 = n2 + 0.07 + offset; break;
    case B: n2 = n2 + 0.13 + offset; break;
  }

  double nn = n1 / n2;
  double ddn = ray.dir.dotProduct(nl);
  double cos = 1 - nn * nn * (1 - ddn*ddn);

  /* Determine transmitted ray direction */
  Vector refraction_dir = (ray.dir * nn - normal * (ddn * nn + sqrt(cos)));

  return radiance(Ray(hitpoint, refraction_dir), wave);
}

int main(int argc, char *argv[]) {
  int width = 1024;
  int height = 768;
  int samples = (argc == 2) ? atoi(argv[1]) : 8;

  double aperture = 2.4;
  double focal_length = 96.5;

  objects.push_back(&light);

  vector<Triangle> layers;

  for (int i = 0; i < 20; i++) {
    const auto h = Vector( 0.0, 250.0,  0.0);
    const auto w = Vector(250.0,  0.0,  0.0);
    const auto d = Vector( 0.0,  0.0, 0.5);

    const auto color = Color(1.0, 1.0, 10.0);

    const auto position = Vector(-100, -100, 100 + i);

    layers.push_back(Triangle(i * d + position, w + h,     h, Color(), color, REFR)); // Front
    layers.push_back(Triangle(i * d + position,     w, w + h, Color(), color, REFR)); // Front
  }

  for (const auto& layer : layers) {
    objects.push_back(&layer);
  }

  /* Set camera origin and viewing direction (negative z direction) */
  Ray camera(Vector(50.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).normalize());
  Vector focal_point = camera.org + camera.dir * focal_length;

  /* Image edge vectors for pixel sampling */
  Vector cx = Vector(width * 0.5135 / height);
  Vector cy = (cx.crossProduct(camera.dir)).normalize() * 0.5135;

  /* Final rendering */
  Image img(width, height);

  /* Loop over image rows */
  for (int y = 0; y < height; y++) {
    cout << "\rRendering (" << samples * 4 << " spp) " << (100.0 * y / (height - 1)) << "%     ";
    srand(y * y * y);

    /* Loop over row pixels */
    #pragma omp parallel for schedule(dynamic, 1)
    for (int x = 0; x < width; x++) {
      img.setColor(x, y, Color());

      /* 2x2 subsampling per pixel */
      for (int sy = 0; sy < 2; sy++) {
        for (int sx = 0; sx < 2; sx++) {
          Color accumulated_radiance = Color();

          /* Compute radiance at subpixel using multiple samples */
          for (int s = 0; s < samples; s++) {
            Vector dir = camera.dir;
            Vector lens_sample_point = Vector(0.0, 0.0, 0.0);

            if (aperture > 0.0) {
              // Generate random sample on circular lens.
              float random_radius = drand48();
              float random_angle = drand48();
              Vector lens_sample_point = aperture * Vector(sqrt(random_radius) * cos(2 * M_PI * random_angle), sqrt(random_radius) * sin(2 * M_PI * random_angle), 0);

              dir = (dir + (focal_point - (camera.org + lens_sample_point)).normalize()).normalize();
            }

            double dx = non_uniform_filter_sample();
            double dy = non_uniform_filter_sample();

            // Ray direction into scene from camera through sample.
            dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
                  cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) +
                  dir;

            auto start = camera.org + dir * 130;

            dir = dir.normalize();

            auto ray = Ray(start + lens_sample_point, dir);


            /* Accumulate radiance */
            accumulated_radiance += Color(
              radiance(ray, R).x,
              radiance(ray, G).y,
              radiance(ray, B).z
            ) / samples;
          }

          accumulated_radiance = accumulated_radiance.clamp(0.0, 1.0) / 4.0;

          img.addColor(x, y, accumulated_radiance);
        }
      }
    }
  }

  cout << endl;
  img.save("image");
}
