/******************************************************************
*
* PathTracing.cpp
*
* Description: This program demonstrates global illumination rendering
* based on the path tracing method. The intergral in the rendering
* equation is approximated via Monte-Carlo integration; explicit
* direct lighting is included to improve quality; the rendered image
* is saved in PPM format.
*
* The code is largely based on the software smallpt by Kevin Beason,
* released under the MIT License.
*
* Advanced Computer Graphics Proseminar WS 2017
*
* Interactive Graphics and Simulation Group
* Department of Computer Science
* University of Innsbruck
*
*******************************************************************/

/* Standard includes */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "Sphere.h"
#include "Triangle.h"
#include "TriangleMeshLoader.h"

#include "../shared/Vector.h"
#include "../shared/Ray.h"
#include "../shared/Image.h"
#include "../shared/macro.h"

using namespace std;

/******************************************************************
* Hard-coded scene definition: the geometry is composed of spheres
* (i.e. Cornell box walls are part of very large spheres).
* These are defined by:
* - radius, center
* - emitted light (light sources), surface reflectivity (~color),
*   material
*******************************************************************/

vector<Sphere> spheres = {
  Sphere(16.5, Vector(27, 16.5, 47), Color(), Color(1.0, 1.0, 1.0),  SPEC), /* Mirror sphere */
  Sphere(16.5, Vector(73, 16.5, 78), Color(), Color(1.0, 1.0, 1.0),  REFR), /* Glas sphere */

  Sphere(1.5, Vector(50, 81.6 - 16.5, 81.6), Color(4, 4, 4) * 100, Color(), DIFF), /* Light */
};

vector<const SceneObject*> objects = vector<const SceneObject*>();

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

Vector randomDirection(Vector direction, double cos_a_max) {
  // Set up local orthogonal coordinate system u, v, w.
  Vector w = direction;
  Vector u = fabs(w.x) > 0.1 ? Vector(0.0, 1.0, 0.0) : Vector(1.0, 0.0, 0.0);
  u = (u.crossProduct(w)).normalize();
  Vector v = w.crossProduct(u);

  double eps1 = drand48();
  double eps2 = drand48();
  double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
  double sin_a = sqrt(1.0 - cos_a * cos_a);
  double phi = 2.0 * M_PI * eps2;

  Vector l = u * cos(phi) * sin_a +
             v * sin(phi) * sin_a +
             w * cos_a;

  return l.normalize();
}

Vector perfectReflection(Vector dir, Vector normal) {
  return dir - normal * 2 * normal.dotProduct(dir);
}

/******************************************************************
* Recursive path tracing for computing radiance via Monte-Carlo
* integration; only considers perfectly diffuse, specular or
* transparent materials;
* after 5 bounces Russian Roulette is used to possibly terminate rays;
* emitted light from light source only included on first direct hit
* (possibly via specular reflection, refraction), controlled by
* parameter E = 0/1;
* on diffuse surfaces light sources are explicitely sampled;
* for transparent objects, Schlick's approximation is employed;
* for first 3 bounces obtain reflected and refracted component,
* afterwards one of the two is chosen randomly
*******************************************************************/
Color radiance(const Ray &ray, int depth, int E, double aperture, double focal_length) {
  depth++;

  double t;
  size_t id = 0;

  if (!intersect(ray, t, id))   /* No intersection with scene */
    return Color(0.0, 0.0, 0.0);

  const SceneObject* obj = objects[id];

  Vector hitpoint = ray.org + ray.dir * t;    /* Intersection point */

  if (depth == 1 && aperture != 0 && focal_length != 0) {
    double dof = 80;

    Vector focal_point = ray.org - Vector(0, 0, focal_length);
    bool behindDepthOfField = hitpoint.z < focal_point.z - dof;
    bool beforeDepthOfField = hitpoint.z > focal_point.z + dof;

    if (beforeDepthOfField || behindDepthOfField) {
      double blur_factor;

      if (beforeDepthOfField) {
        blur_factor = clamp(hitpoint.z - (focal_point.z + dof), 0.0, dof) / dof;
      } else {
        blur_factor = clamp(hitpoint.z - (focal_point.z - dof), -dof, 0.0) / -dof;
      }

      double cos_a_max = cos(0.005 + blur_factor * 0.00002);
      Vector l = randomDirection(ray.dir, cos_a_max);

      return radiance(Ray(ray.org, l), 0, E, 0, 0);
    }

    // // https://en.wikipedia.org/wiki/Circle_of_confusion
    // double object_distance = fabs(hitpoint.z - ray.org.z); // Object Distance
    // double image_distance = focal_length * object_distance / (object_distance - focal_length);
    // double focused_object_distance = focal_length * image_distance / (image_distance - focal_length);
    //
    // double m = image_distance / focused_object_distance;
    // double C = aperture * fabs(object_distance - focused_object_distance) / object_distance;
    // double c = C * m;
    //
    // double N = focal_length / aperture; // https://en.wikipedia.org/wiki/F-number
    //
    // double dof = 2 * N * c * (m + 1) / (pow(m, 2) - pow(N * c / focal_length, 2));
    //
    // double cos_a_max = cos(atan(c / 2 / image_distance));
    // Vector l = randomDirection(ray.dir, cos_a_max);
    //
    // return radiance(Ray(ray.org, l), 0, E, 0, 0);
  }

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

  Color col = obj->color;

  /* Maximum RGB reflectivity for Russian Roulette */
  double p = clamp(col.max(), 0.0, 0.999);

  /* After 5 bounces or if max reflectivity is zero */
  if (depth > 5 || p == 0) {
    /* Russian Roulette */
    if (drand48() >= p)
      return obj->emission * E; /* No further bounces, only return potential emission */

    col = col * (1.0 / p); /* Scale estimator to remain unbiased */
  }

  if (obj->refl == DIFF) {
    /* Compute random reflection vector on hemisphere */
    double r1 = 2.0 * M_PI * drand48();
    double r2 = drand48();
    double r2s = sqrt(r2);

    /* Set up local orthogonal coordinate system u,v,w on surface */
    Vector w = nl;
    Vector u = fabs(w.x) > 0.1 ? Vector(0.0, 1.0, 0.0) : Vector(1.0, 0.0, 0.0);
    u = u.crossProduct(w).normalize();

    Vector v = w.crossProduct(u);

    /* Random reflection vector d */
    Vector d = (u * cos(r1) * r2s +
                v * sin(r1) * r2s +
                w * sqrt(1 - r2)).normalize();

    /* Explicit computation of direct lighting */
    Vector e;
    for (size_t i = 0; i < objects.size(); i++) {
      const SceneObject* lightSource = objects[i];
      if (lightSource->emission.x <= 0 && lightSource->emission.y <= 0 && lightSource->emission.z <= 0)
          continue; /* Skip objects that are not light sources */

      if (!lightSource->isSphere) {
        cerr << "Warning: Only spherical light sources are implemented." << endl;
        continue;
      }

      const Sphere* sphere = static_cast<const Sphere*>(lightSource);

      /* Randomly sample spherical light source from surface intersection */

      // Create random sample direction towards spherical light source.
      double cos_a_max = sqrt(1.0 - pow(sphere->radius, 2) /
                              (hitpoint - sphere->position).dotProduct(hitpoint - sphere->position));
      Vector l = randomDirection(sphere->position - hitpoint, cos_a_max);

      /* Shoot shadow ray, check if intersection is with light source */
      if (intersect(Ray(hitpoint, l), t, id) && id == i) {
        double omega = 2 * M_PI * (1 - cos_a_max);

        /* Add diffusely reflected light from light source; note constant BRDF 1 / π */
        e += col.entrywiseProduct(sphere->emission * l.dotProduct(nl) * omega) / M_PI;
      }
    }

    /* Return potential light emission, direct lighting, and indirect lighting (via
        recursive call for Monte-Carlo integration */
    return obj->emission * E + e + col.entrywiseProduct(radiance(Ray(hitpoint, d), depth, 0, aperture, focal_length));
  } else if (obj->refl == GLOS) {
    double cos_a_max = cos(0.10);
    Vector l = randomDirection(nl, cos_a_max);

    return obj->emission + col.entrywiseProduct(radiance(Ray(hitpoint, l), depth, 1, aperture, focal_length));
  }

  if (obj->refl == SPEC) {
    // Return light emission mirror reflection (via recursive call using perfect reflection vector).
    Ray reflection_ray = Ray(hitpoint, perfectReflection(ray.dir, normal));
    return obj->emission + col.entrywiseProduct(radiance(reflection_ray, depth, 1, aperture, focal_length));
  }

  Vector ray_dir = ray.dir;

  /* Otherwise object transparent, i.e. assumed dielectric glass material */
  double nc = 1;                        /* Index of refraction of air (approximately) */
  double nt = 1.5;                      /* Index of refraction of glass (approximately) */

  if (obj->refl == TRAN) {
    double cos_a_max = cos(0.2);
    ray_dir = randomDirection(ray.dir, cos_a_max);
    nt = 1.15;
  }

  bool into = normal.dotProduct(nl) > 0;       /* Bool for checking if ray from outside going in */
  double nnt = into ? nc / nt : nt / nc; /* Set ratio depending on hit from inside or outside */

  double ddn = ray_dir.dotProduct(nl);
  double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

  Ray reflection_ray = Ray(hitpoint, perfectReflection(ray_dir, normal)); // Perfect reflection.

  // Check for total internal reflection, if so only reflect.
  if (cos2t < 0)
    return obj->emission + col.entrywiseProduct(radiance(reflection_ray, depth, 1, aperture, focal_length));

  /* Otherwise reflection and/or refraction occurs */
  Vector transmission_direction;

  // Determine transmitted ray direction for refraction.
  if(into) {
    transmission_direction = ray_dir * nnt - normal * (ddn * nnt + sqrt(cos2t));
  } else {
    transmission_direction = ray_dir * nnt + normal * (ddn * nnt + sqrt(cos2t));
  }

  transmission_direction = transmission_direction.normalize();

  /* Cosine of correct angle depending on outside/inside */
  double c = into ? 1 + ddn : 1 - transmission_direction.dotProduct(normal);

  auto schlickReflectance = [](double n1, double n2, double c) {
    const double R0 = pow((n1 - n2) / (n1 + n2), 2);
    return R0 + (1 - R0) * pow(c, 5);
  };

  /* Compute Schlick's approximation of Fresnel equation */
  double Re = schlickReflectance(nc, nt, c);   // Reflectance
  double Tr = 1 - Re;                     // Transmittance

  Ray transmission_ray = Ray(hitpoint, transmission_direction);

  // Initially, use both reflection and trasmission.
  if (depth < 3)
    return obj->emission + col.entrywiseProduct(radiance(reflection_ray, depth, 1, aperture, focal_length) * Re +
                                                radiance(transmission_ray, depth, 1, aperture, focal_length) * Tr);

  // Probability for selecting reflectance or transmittance */
  double P = 0.25 + 0.5 * Re;
  double RP = Re / P;         /* Scaling factors for unbiased estimator */
  double TP = Tr / (1 - P);

  /* Russian Roulette */
  if (drand48() < P)
    return obj->emission + col.entrywiseProduct(radiance(reflection_ray, depth, 1, aperture, focal_length) * RP);

  return obj->emission + col.entrywiseProduct(radiance(transmission_ray, depth, 1, aperture, focal_length) * TP);
}


/******************************************************************
* Main routine: Computation of path tracing image (2x2 subpixels)
* Key parameters
* - Image dimensions: width, height
* - Number of samples per subpixel (non-uniform filtering): samples
* Rendered result saved as PPM image file
*******************************************************************/

int main(int argc, char *argv[]) {
  int width = 1024;
  int height = 768;
  int samples = (argc == 2) ? atoi(argv[1]) : 12;

  double aperture = 1.4;
  double focal_length = 40;

  const auto& walls = TriangleMeshLoader::loadTriangleMesh("meshes/Walls.trim");

  if(walls) {
    for (auto &wall : *walls) {
      objects.push_back(&wall);
    }
  }

  const auto& cuboid = TriangleMeshLoader::loadTriangleMesh("meshes/Cuboid.trim");

  if(cuboid) {
    for (auto &cb : *cuboid) {
      objects.push_back(&cb);
    }
  }

  for (auto &sphere : spheres) {
    objects.push_back(&sphere);
  }

  /* Set camera origin and viewing direction (negative z direction) */
  Ray camera(Vector(50.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).normalize());

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
    #pragma omp parallel for
    for (int x = 0; x < width; x++) {
      img.setColor(x, y, Color());

      /* 2x2 subsampling per pixel */
      for (int sy = 0; sy < 2; sy++) {
        for (int sx = 0; sx < 2; sx++) {
          Color accumulated_radiance = Color();

          /* Compute radiance at subpixel using multiple samples */
          for (int s = 0; s < samples; s++) {
            double dx = non_uniform_filter_sample();
            double dy = non_uniform_filter_sample();

            /* Ray direction into scene from camera through sample */
            Vector dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
                         cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) +
                         camera.dir;

            /* Extend camera ray to start inside box */
            Vector start = camera.org + dir * 130.0;

            dir = dir.normalize();

            /* Accumulate radiance */
            accumulated_radiance = accumulated_radiance +
              radiance(Ray(start, dir), 0, 1, aperture, focal_length) / samples;
          }

          accumulated_radiance = accumulated_radiance.clamp() * 0.25;

          img.addColor(x, y, accumulated_radiance);
        }
      }
    }
  }

  cout << endl;
  img.save("image");
}
