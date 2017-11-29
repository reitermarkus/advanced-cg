/******************************************************************
 *
 * Radiosity.cpp
 *
 * Description: This file demonstrates global illumination rendering
 * based on the radiosity method. The geometry is divided into patches
 * for which the form factors are determined employing Monte Carlo
 * integration. Radiosity values for the patches are computed with
 * an iterative solver. The final image (i.e. radiance) is obtained
 * via tracing rays into the scene. Two output files are saved -
 * one with constant shading of patches and one with bicubic color
 * interpolation.
 *
 * The code is extended from software by user Hole and Kevin Beason,
 * released under the MIT License.
 *
 * http://kagamin.net/hole/license.txt
 * http://kagamin.net/hole/smallpt-license.txt
 *
 * Advanced Computer Graphics Proseminar WS 2015
 *
 * Interactive Graphics and Simulation Group
 * Institute of Computer Science
 * University of Innsbruck
 *
 *******************************************************************/

/* Standard includes */
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <chrono>

#include "Vector.h"
#include "Image.h"
#include "Triangle.h"
#include "Ray.h"
#include "PatchTriangle.h"

#include "macro.h"

#undef M_PI
const double M_PI = atan(1) * 4;

using namespace std;
using namespace std::chrono;

static map<Triangle*, vector<map<Triangle*, vector<double>>>> form_factor;
static int patch_num = 0;


const Color backgroundColor(0.0, 0.0, 0.0);

/******************************************************************
 * Hard-coded scene definition: the geometry is composed of triangles.
 * These are defined by:
 * - vector to corner(origin), b_rel, c_rel
 * - emitted light energy (light sources), surface reflectivity (~color)
 *******************************************************************/
vector<Triangle> tris = {
  /* Cornell Box walls */
  Triangle(Vector(  0.0,  0.0,   0.0), Vector( 100.0, 0.0,    0.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.75, 0.75, 0.75)), // Back:   bottom-left
  Triangle(Vector(100.0, 80.0,   0.0), Vector(-100.0, 0.0,    0.0), Vector(0.0, -80.0,    0.0), Color(), Color(0.75, 0.75, 0.75)), // Back:   top-right
  Triangle(Vector(  0.0,  0.0, 170.0), Vector( 100.0, 0.0,    0.0), Vector(0.0,   0.0, -170.0), Color(), Color(0.75, 0.75, 0.75)), // Bottom: front-left
  Triangle(Vector(100.0,  0.0,   0.0), Vector(-100.0, 0.0,    0.0), Vector(0.0,   0.0,  170.0), Color(), Color(0.75, 0.75, 0.75)), // Bottom: back-right
  Triangle(Vector(  0.0, 80.0,   0.0), Vector( 100.0, 0.0,    0.0), Vector(0.0,   0.0,  170.0), Color(), Color(0.75, 0.75, 0.75)), // Top:    back-left
  Triangle(Vector(100.0, 80.0, 170.0), Vector(-100.0, 0.0,    0.0), Vector(0.0,   0.0, -170.0), Color(), Color(0.75, 0.75, 0.75)), // Top:    front-right
  Triangle(Vector(  0.0,  0.0, 170.0), Vector(   0.0, 0.0, -170.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.75, 0.25, 0.25)), // Left:   front-bottom
  Triangle(Vector(  0.0, 80.0,   0.0), Vector(   0.0, 0.0,  170.0), Vector(0.0, -80.0,    0.0), Color(), Color(0.75, 0.25, 0.25)), // Left:   back-top
  Triangle(Vector(100.0,  0.0,   0.0), Vector(   0.0, 0.0,  170.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.25, 0.25, 0.75)), // Right:  back-bottom
  Triangle(Vector(100.0, 80.0, 170.0), Vector(   0.0, 0.0, -170.0), Vector(0.0, -80.0,    0.0), Color(), Color(0.25, 0.25, 0.75)), // Right:  front-top
  Triangle(Vector(100.0,  0.0, 170.0), Vector(-100.0, 0.0,    0.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.25, 0.75, 0.25)),  // Front:  bottom-right (not visible)
  Triangle(Vector(  0.0, 80.0, 170.0), Vector( 100.0, 0.0,    0.0), Vector(0.0, -80.0,    0.0), Color(), Color(0.25, 0.75, 0.25)),  // Front:  top-left (not visible)

  /* Area light source on top */
  Triangle(Vector(40.0, 79.99, 65.0), Vector( 20.0, 0.0, 0.0), Vector(0.0, 0.0,  20.0), Color(12, 12, 12), Color(0.75, 0.75, 0.75)), // back-left
  Triangle(Vector(60.0, 79.99, 85.0), Vector(-20.0, 0.0, 0.0), Vector(0.0, 0.0, -20.0), Color(12, 12, 12), Color(0.75, 0.75, 0.75)), // front-right

  /* Cuboid in room */
  Triangle(Vector(30.0,  0.0, 100.0), Vector(  0.0, 0.0, -20.0), Vector(0.0,  40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Right: front-bottom
  Triangle(Vector(30.0, 40.0,  80.0), Vector(  0.0, 0.0,  20.0), Vector(0.0, -40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Right: back-top
  Triangle(Vector(10.0,  0.0,  80.0), Vector(  0.0, 0.0,  20.0), Vector(0.0,  40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Left:  back-bottom
  Triangle(Vector(10.0, 40.0, 100.0), Vector(  0.0, 0.0, -20.0), Vector(0.0, -40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Left:  front-top
  Triangle(Vector(10.0,  0.0, 100.0), Vector( 20.0, 0.0,   0.0), Vector(0.0,  40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Front: bottom-left
  Triangle(Vector(30.0, 40.0, 100.0), Vector(-20.0, 0.0,   0.0), Vector(0.0, -40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Front: top-right
  Triangle(Vector(30.0,  0.0,  80.0), Vector(-20.0, 0.0,   0.0), Vector(0.0,  40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Back:  bottom-right
  Triangle(Vector(10.0, 40.0,  80.0), Vector( 20.0, 0.0,   0.0), Vector(0.0, -40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Back:  top-left
  Triangle(Vector(10.0, 40.0, 100.0), Vector( 20.0, 0.0,   0.0), Vector(0.0,   0.0, -20.0), Color(), Color(0.75, 0.75, 0.75)), // Top:   front-left
  Triangle(Vector(30.0, 40.0,  80.0), Vector(-20.0, 0.0,   0.0), Vector(0.0,   0.0,  20.0), Color(), Color(0.75, 0.75, 0.75)), // Top:   back-right
};

/******************************************************************
 * Check for closest intersection of a ray with the scene;
 * Returns true if intersection is found, as well as ray parameter
 * of intersection and id of intersected object
 *******************************************************************/
bool intersectScene(const Ray &ray, double *t, int *id, Vector *normal) {
  const int n = tris.size();
  *t = 1e20;
  *id = -1;

  for (int i = 0; i < n; i++) {
    double d = tris[i].intersect(ray);
    if (d > 0.0 && d < *t) {
      *t = d;
      *id = i;
      *normal = tris[i].normal;
    }
  }
  return *t < 1e20;
}

/******************************************************************
 * Determine all form factors for all pairs of patches (of all
 * rectangles);
 * Evaluation of integrals in form factor equation is done via
 * Monte Carlo integration; samples are uniformly distributed and
 * equally weighted;
 * Computation accelerated by exploiting symmetries of form factor
 * estimation;
 *******************************************************************/
void calculateFormFactors(const int divisions, const int mc_sample) {
  /* Total number of patches in scene */
  const int n = tris.size();
  for (auto &tri : tris) {
    tri.init_patches(divisions);
    patch_num += tri.patch.size();
  }
  unsigned long form_factor_num = pow(patch_num, 2);

  cout << "Number of triangles: " << n << endl;
  cout << "Number of patches: " << patch_num << endl;
  cout << "Number of form factors: " << form_factor_num << endl;

  for (auto &tri_a : tris) {
    form_factor[&tri_a] = vector<map<Triangle*, vector<double>>>(tri_a.patch.size());
    for (size_t p = 0; p < tri_a.patch.size(); p++) {
      form_factor[&tri_a][p] = map<Triangle*, vector<double>>();
      for (auto &tri_b : tris) {
        form_factor[&tri_a][p][&tri_b] = vector<double>(tri_b.patch.size());
        fill(form_factor[&tri_a][p][&tri_b].begin(), form_factor[&tri_a][p][&tri_b].end(), 0.0);
      }
    }
  }

  map<Triangle*, vector<double>> patch_area;

  /* Precompute patch areas, assuming same size for each triangle */
  for (auto &tri : tris) {
    patch_area[&tri] = vector<double>(tri.patch.size());
    const auto area = tri.area / tri.patch.size();
    fill(patch_area[&tri].begin(), patch_area[&tri].end(), area);
  }

  /* Loop over all triangles in scene */
  for (int i = 0; i < n; i++) {
    cout << i << " ";

    /* Loop over all patches in rectangle i */
    #pragma omp parallel for
    for (unsigned long p_i = 0; p_i < tris[i].patch.size(); p_i++) {
      if (p_i % tris[i].divisions == 0) {
        cout << "*" << flush;
      }

      /* Loop over all triangles in scene for triangles i */
      for (int j = i + 1; j < n; j++) {
        /* Loop over all patches in rectangle j */
        for (unsigned long p_j = 0; p_j < tris[j].patch.size(); p_j++) {
          double F = 0;

          /* Monte Carlo integration of form factor double integral */

          /* Uniform PDF for Monte Carlo (1/Ai)x(1/Aj) */
          const double pdf =
              (1.0 / patch_area[&tris[i]][p_i]) *
              (1.0 / patch_area[&tris[j]][p_j]);

          for (auto ii = 0; ii < mc_sample; ii++) {
            PatchTriangle t_i = tris[i].subTriangles[p_i];
            PatchTriangle t_j = tris[j].subTriangles[p_j];

            const Vector xi = Triangle::sample(t_i.a, t_i.b, t_i.c);
            const Vector xj = Triangle::sample(t_j.a, t_j.b, t_j.c);

            /* Check for visibility between sample points */
            const Vector ij = (xj - xi).normalize();

            double t;
            int id;
            Vector normal;
            if (intersectScene(Ray(xi, ij), &t, &id, &normal) && id != j) {
              continue; /* If intersection with other triangle */
            }

            /* Cosines of angles beteen normals and ray inbetween */
            const double d0 = tris[i].normal.dotProduct(ij);
            const double d1 = tris[j].normal.dotProduct(-1.0 * ij);

            /* Continue if patches facing each other */
            if (d0 > 0.0 && d1 > 0.0) {
              /* Sample form factor */
              const double K = d0 * d1 / (M_PI * (xj - xi).lengthSquared());

              /* Add weighted sample to estimate */
              F += K / pdf;
            }
          }

          /* Divide by number of samples */
          F /= mc_sample;

          form_factor[&tris[i]][p_i][&tris[j]][p_j] = F;
          form_factor[&tris[j]][p_j][&tris[i]][p_i] = F;
        }
      }
    }

    cout << endl;
  }

  /* Divide by area to get final form factors */
  for (auto &tri_a : tris) {
    for (size_t p_a = 0; p_a < tri_a.patch.size(); p_a++) {
      for (auto &tri_b : tris) {
        for (size_t p_b = 0; p_b < tri_b.patch.size(); p_b++) {
          const auto area = patch_area[&tri_a][p_a];
          form_factor[&tri_a][p_a][&tri_b][p_b] = clamp(form_factor[&tri_a][p_a][&tri_b][p_b] / area, 0.0, 1.0);
        }
      }
    }
  }
}

unordered_map<Vector, unordered_map<Vector, Color>> calculateVertexColors() {
  auto vertex_colors = unordered_map<Vector, unordered_map<Vector, Color>>();
  auto vertex_counts = unordered_map<Vector, unordered_map<Vector, int>>();

  for (auto &tri : tris) {
    vertex_colors.emplace(tri.normal, unordered_map<Vector, Color>());
    vertex_counts.emplace(tri.normal, unordered_map<Vector, int>());
  }

  for (auto &tri : tris) {
    for (size_t p = 0; p < tri.patch.size(); p++) {
      vertex_colors.at(tri.normal).emplace(tri.subTriangles[p].a, Color(0, 0, 0));
      vertex_colors.at(tri.normal).emplace(tri.subTriangles[p].b, Color(0, 0, 0));
      vertex_colors.at(tri.normal).emplace(tri.subTriangles[p].c, Color(0, 0, 0));
      vertex_counts.at(tri.normal).emplace(tri.subTriangles[p].a, 0);
      vertex_counts.at(tri.normal).emplace(tri.subTriangles[p].b, 0);
      vertex_counts.at(tri.normal).emplace(tri.subTriangles[p].c, 0);
    }
  }

  for (auto &tri : tris) {
    for (size_t p = 0; p < tri.patch.size(); p++) {
      vertex_counts.at(tri.normal).at(tri.subTriangles[p].a) += 1;
      vertex_counts.at(tri.normal).at(tri.subTriangles[p].b) += 1;
      vertex_counts.at(tri.normal).at(tri.subTriangles[p].c) += 1;
    }
  }

  for (auto &tri : tris) {
    for (size_t p = 0; p < tri.patch.size(); p++) {
      vertex_colors.at(tri.normal).at(tri.subTriangles[p].a) += tri.patch[p];
      vertex_colors.at(tri.normal).at(tri.subTriangles[p].b) += tri.patch[p];
      vertex_colors.at(tri.normal).at(tri.subTriangles[p].c) += tri.patch[p];
    }
  }

  for (auto &entry_1 : vertex_colors) {
    auto normal = entry_1.first;
    auto vertices = entry_1.second;

    for (auto &entry_2 : vertices) {
      auto vertex = entry_2.first;
      vertex_colors.at(normal).at(vertex) /= vertex_counts.at(normal).at(vertex);
    }
  }

  return vertex_colors;
}

/******************************************************************
 * Iterative computation of radiosity via Gathering; i.e. solution
 * using Gauss-Seidel iteration - reuse already computed values;
 * run-time O(n^2)
 *******************************************************************/

void calculateRadiosity() {
  for (auto &tri_a : tris) {
    #pragma omp parallel for
    for (unsigned long p_a = 0; p_a < tri_a.patch.size(); p_a++) {
      Color B;

      for (auto &tri_b : tris) {
        for (unsigned long p_b = 0; p_b < tri_b.patch.size(); p_b++) {
          const double Fij = form_factor[&tri_a][p_a][&tri_b][p_b];

          /* Add form factor multiplied with radiosity of previous step */
          if (Fij > 0.0)
            B = B + Fij * tri_b.patch[p_b];
        }
      }

      /* Multiply sum with color of patch and add emission */
      B = tri_a.color.entrywiseProduct(B) + tri_a.emission;

      /* Store overall patch radiosity of current iteration */
      tri_a.patch[p_a] = B;
    }
  }
}

/******************************************************************
 * Compute radiance from radiosity by shooting rays into the scene;
 * Radiance directly proportional to radiosity for assumed diffuse
 * emitters/surfaces (multiply by PI);
 * At intersections either constant patch color is returned or a
 * smoothly interpolated color of 4x4 neighboring patches
 *******************************************************************/

static int findPatchIndex(const Triangle &triangle, const Ray &ray) {
  for (unsigned long p = 0; p < triangle.subTriangles.size(); p++) {
    if (triangle.subTriangles[p].intersect(ray)) {
      return p;
    }
  }

  return -1;
}

pair<Color, Color> radiance(const Ray &ray, unordered_map<Vector, unordered_map<Vector, Color>> &vertex_colors) {
  double t;
  int id;
  Vector normal;

  /* Find intersected triangle. */
  if (!intersectScene(ray, &t, &id, &normal)) {
    return make_pair(backgroundColor, backgroundColor);
  }

  /* Find intersected patch. */
  const Triangle &obj = tris[id];
  int idx = findPatchIndex(obj, ray);

  const Vector hitpoint = ray.org + t * ray.dir;

  Vector bary = obj.subTriangles[idx].barycentricCoordinatesAt(hitpoint);

  Color a = vertex_colors[obj.normal][obj.subTriangles[idx].a];
  Color b = vertex_colors[obj.normal][obj.subTriangles[idx].b];
  Color c = vertex_colors[obj.normal][obj.subTriangles[idx].c];

  Color interpolated = a * bary.x + b * bary.z + c * bary.y;

  return make_pair(obj.patch[idx], interpolated);
}

/******************************************************************
 * Main routine: Computation of radiosity image
 * Key parameters
 * - Image dimensions: width, height
 * - Number of samples for antialiasing (non-uniform filter): samples
 * - Number of patches along edges a,b: patches_a, patches_b
 * - Number of uniform samples per patch edge: MC_samples
 * - Number of iterations for iterative solver: iterations
 * Rendered result saved as PPM image file
 *******************************************************************/

int main(void) {
  duration<double, milli> elapsed;
  auto total_bench = high_resolution_clock::now();

  int divisions = 10;
  int mc_samples = 10;

  cout << "Calculating form factors" << endl;
  auto form_factor_bench = high_resolution_clock::now();
  calculateFormFactors(divisions, mc_samples);

  elapsed = high_resolution_clock::now() - form_factor_bench;
  auto form_factor_bench_elapsed = elapsed.count();

  /* Iterative solution of radiosity linear system */
  cout << "Calculating radiosity ..." << endl;
  int iterations = 40;
  auto radiosity_bench = high_resolution_clock::now();
  for (int i = 0; i < iterations; i++) {
    cout << i << " ";
    calculateRadiosity();
  }
  cout << endl;

  elapsed = high_resolution_clock::now() - radiosity_bench;
  auto radiosity_bench_elapsed = elapsed.count();


  cout << "Calculating vertex colors ..." << endl;
  auto vertex_colors = calculateVertexColors();


  cout << "Rendering images ..." << endl;

  int width = 640;
  int height = 480;

  int samples = 4;

  /* Set camera origin and viewing direction (negative z direction) */
  Ray camera(Vector(50.0, 52.0, 295.6),
             Vector(0.0, -0.042612, -1.0).normalize());

  /* Image edge vectors for pixel sampling */
  Vector cx = Vector(width * 0.5135 / height);
  Vector cy = (cx.crossProduct(camera.dir)).normalize() * 0.5135;

  /* Two final renderings; one with constant, one with interpolated patch colors
   */
  Image img(width, height);
  Image img_interpolated(width, height);

  auto rendering_bench = high_resolution_clock::now();

  /* Loop over image rows */
  for (int y = 0; y < height; y++) {
    cout << "\rRendering (" << samples * 4 << " spp) "
         << (100.0 * y / (height - 1)) << "%     ";
    srand(y * y * y);

    /* Loop over row pixels */
    #pragma omp parallel for
    for (int x = 0; x < width; x++) {
      img.setColor(x, y, Color());
      img_interpolated.setColor(x, y, Color());

      /* 2x2 subsampling per pixel */
      for (int sy = 0; sy < 2; sy++) {
        for (int sx = 0; sx < 2; sx++) {
          Color accumulated_radiance = Color();
          Color accumulated_radiance2 = Color();

          /* Computes radiance at subpixel using multiple samples */
          for (int s = 0; s < samples; s++) {
            auto nu_filter_samples = [] {
              /* Transform uniform into non-uniform filter samples */
              auto r = 2.0 * drand48();
              return r < 1.0 ? (sqrt(r) - 1.0) : (1.0 - sqrt(2.0 - r));
            };

            double dx = nu_filter_samples();
            double dy = nu_filter_samples();

            /* Ray direction into scene from camera through sample */
            Vector dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
                         cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) +
                         camera.dir;

            /* Extend camera ray to start inside box */
            Vector start = camera.org + dir * 130.0;

            auto colors = radiance(Ray(start, dir.normalize()), vertex_colors);

            /* Determine constant radiance */
            accumulated_radiance =
                accumulated_radiance +
                colors.first / samples;

            /* Determine interpolated radiance */
            accumulated_radiance2 =
                accumulated_radiance2 +
                colors.second / samples;
          }

          img.addColor(x, y, accumulated_radiance);
          img_interpolated.addColor(x, y, accumulated_radiance2);
        }
      }
    }
  }

  cout << endl;

  elapsed = high_resolution_clock::now() - rendering_bench;
  auto rendering_bench_elapsed = elapsed.count();


  cout << "Saving images ..." << endl;
  auto image_bench = high_resolution_clock::now();

  img.save(string("image_patches.ppm"));
  img_interpolated.save(string("image_smooth.ppm"));

  elapsed = high_resolution_clock::now() - image_bench;
  auto image_bench_elapsed = elapsed.count();

  elapsed = high_resolution_clock::now() - total_bench;
  auto total_bench_elapsed = elapsed.count();

  printf("┌───────────────┬───────────┐\n");
  printf("│ Form Factors  │ %6ld ms │\n", (unsigned long)form_factor_bench_elapsed);
  printf("├───────────────┼───────────┤\n");
  printf("│ Radiosity     │ %6ld ms │\n", (unsigned long)radiosity_bench_elapsed);
  printf("├───────────────┼───────────┤\n");
  printf("│ Rendering     │ %6ld ms │\n", (unsigned long)rendering_bench_elapsed);
  printf("├───────────────┼───────────┤\n");
  printf("│ Saving Images │ %6ld ms │\n", (unsigned long)image_bench_elapsed);
  printf("┢━━━━━━━━━━━━━━━╈━━━━━━━━━━━┪\n");
  printf("┃ Total         ┃ %6ld ms ┃\n", (unsigned long)total_bench_elapsed);
  printf("┗━━━━━━━━━━━━━━━┻━━━━━━━━━━━┛\n");

}
