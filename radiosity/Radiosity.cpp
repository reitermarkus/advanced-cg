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
#include <vector>

#include "Vector.h"
#include "ColorUtils.h"
#include "Image.h"
#include "Rectangle.h"
#include "Triangle.h"
#include "Ray.h"

#undef M_PI
const double M_PI = atan(1) * 4;

using namespace std;

static map<Rectangle*, vector<map<Rectangle*, vector<double>>>> form_factor;
static int patch_num = 0;

const Color BackgroundColor(0.0, 0.0, 0.0);

/******************************************************************
 * Hard-coded scene definition: the geometry is composed of triangles.
 * These are defined by:
 * - vector to corner(origin), edge a, edge b
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
  Triangle(Vector(100.0,  0.0, 170.0), Vector(-100.0, 0.0,    0.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.0,  1.0,  0.0)),  // Front:  bottom-right (not visible)
  Triangle(Vector(  0.0, 80.0, 170.0), Vector( 100.0, 0.0,    0.0), Vector(0.0, -80.0,    0.0), Color(), Color(0.0,  1.0,  0.0)),  // Front:  top-left (not visible)

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
  Triangle(Vector(10.0,  4.0,  80.0), Vector( 20.0, 0.0,   0.0), Vector(0.0, -40.0,   0.0), Color(), Color(0.75, 0.75, 0.75)), // Back:  top-left
  Triangle(Vector(10.0, 40.0, 100.0), Vector( 20.0, 0.0,   0.0), Vector(0.0,   0.0, -20.0), Color(), Color(0.75, 0.75, 0.75)), // Top:   front-left
  Triangle(Vector(30.0, 40.0,  80.0), Vector(-20.0, 0.0,   0.0), Vector(0.0,   0.0,  20.0), Color(), Color(0.75, 0.75, 0.75)), // Top:   back-right
};

/******************************************************************
 * Hard-coded scene definition: the geometry is composed of rectangles.
 * These are defined by:
 * - vector to corner(origin), edge a, edge b
 * - emitted light energy (light sources), surface reflectivity (~color)
 *******************************************************************/
vector<Rectangle> recs = {
  /* Cornell Box walls */
  Rectangle(Vector(  0.0,  0.0,   0.0), Vector( 100.0, 0.0,    0.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.75, 0.75, 0.75)), /* Back */
  Rectangle(Vector(  0.0,  0.0, 170.0), Vector( 100.0, 0.0,    0.0), Vector(0.0,   0.0, -170.0), Color(), Color(0.75, 0.75, 0.75)), /* Bottom */
  Rectangle(Vector(  0.0, 80.0,   0.0), Vector( 100.0, 0.0,    0.0), Vector(0.0,   0.0,  170.0), Color(), Color(0.75, 0.75, 0.75)), /* Top */
  Rectangle(Vector(  0.0,  0.0, 170.0), Vector(   0.0, 0.0, -170.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.75, 0.25, 0.25)), /* Left */
  Rectangle(Vector(100.0,  0.0,   0.0), Vector(   0.0, 0.0,  170.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.25, 0.25, 0.75)), /* Right */
  Rectangle(Vector(100.0,  0.0, 170.0), Vector(-100.0, 0.0,    0.0), Vector(0.0,  80.0,    0.0), Color(), Color(0.0,  1.0,  0.0)),  /* Front (not visible) */

  /* Area light source on top */
  Rectangle(Vector(40.0, 79.99, 65.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, 20.0), Color(12, 12, 12), Color(0.75, 0.75, 0.75)),

  /* Cuboid in room */
  Rectangle(Vector(30.0,  0.0, 100.0), Vector(  0.0, 0.0, -20.0), Vector(0.0,  40.0,  0.0), Color(), Color(0.75, 0.75, 0.75)), /* Right */
  Rectangle(Vector(10.0,  0.0,  80.0), Vector(  0.0, 0.0,  20.0), Vector(0.0,  40.0,  0.0), Color(), Color(0.75, 0.75, 0.75)), /* Left */
  Rectangle(Vector(10.0,  0.0, 100.0), Vector( 20.0, 0.0,   0.0), Vector(0.0,  40.0,  0.0), Color(), Color(0.75, 0.75, 0.75)), /* Front */
  Rectangle(Vector(30.0,  0.0,  80.0), Vector(-20.0, 0.0,   0.0), Vector(0.0,  40.0,  0.0), Color(), Color(0.75, 0.75, 0.75)), /* Back */
  Rectangle(Vector(10.0, 40.0, 100.0), Vector( 20.0, 0.0,   0.0), Vector(0.0,  0.0, -20.0), Color(), Color(0.75, 0.75, 0.75)), /* Top */
};

/******************************************************************
 * Check for closest intersection of a ray with the scene;
 * Returns true if intersection is found, as well as ray parameter
 * of intersection and id of intersected object
 *******************************************************************/
bool intersectScene(const Ray &ray, double *t, int *id, Vector *normal) {
  const int n = recs.size();
  *t = 1e20;
  *id = -1;

  for (int i = 0; i < n; i++) {
    double d = recs[i].intersect(ray);
    if (d > 0.0 && d < *t) {
      *t = d;
      *id = i;
      *normal = recs[i].normal;
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
void calculateFormFactors(const int a_div_num, const int b_div_num,
                            const int mc_sample) {
  /* Total number of patches in scene */
  const int n = recs.size();
  for (auto &rec : recs) {
    rec.init_patches(a_div_num, b_div_num);
    patch_num += rec.patch.size();
  }
  int form_factor_num = pow(patch_num, 2);

  cout << "Number of rectangles: " << n << endl;
  cout << "Number of patches: " << patch_num << endl;
  cout << "Number of form factors: " << form_factor_num << endl;

  for (auto &rec_a : recs) {
    form_factor[&rec_a] = vector<map<Rectangle*, vector<double>>>(rec_a.patch.size());
    for (size_t p = 0; p < rec_a.patch.size(); p++) {
      form_factor[&rec_a][p] = map<Rectangle*, vector<double>>();
      for (auto &rec_b : recs) {
        form_factor[&rec_a][p][&rec_b] = vector<double>(rec_b.patch.size());
        fill(form_factor[&rec_a][p][&rec_b].begin(), form_factor[&rec_a][p][&rec_b].end(), 0.0);
      }
    }
  }

  map<Rectangle*, vector<double>> patch_area;

  /* Precompute patch areas, assuming same size for each rectangle */
  for (auto &rec : recs) {
    patch_area[&rec] = vector<double>(rec.patch.size());
    const auto area = rec.area / rec.patch.size();
    fill(patch_area[&rec].begin(), patch_area[&rec].end(), area);
  }

  /* Loop over all rectangles in scene */
  for (int i = 0; i < n; i++) {
    cout << i << " ";

    /* Loop over all patches in rectangle i */
    #pragma omp parallel for
    for (int ia = 0; ia < recs[i].a_num; ia++) {
      cout << "*" << flush;
      for (int ib = 0; ib < recs[i].b_num; ib++) {
        const Vector normal_i = recs[i].normal;

        /* Loop over all rectangles in scene for rectangle i */
        for (int j = 0; j < n; j++) {
          const Vector normal_j = recs[j].normal;

          /* Loop over all patches in rectangle j */
          for (int ja = 0; ja < recs[j].a_num; ja++) {
            for (int jb = 0; jb < recs[j].b_num; jb++) {
              /* Do not compute form factors for patches on same rectangle;
                 also exploit symmetry to reduce computation;
                 intemediate values; will be divided by patch area below */
              if (i < j) {
                double F = 0;

                /* Monte Carlo integration of form factor double integral */
                const int Ni = mc_sample, Nj = mc_sample;

                /* Uniform PDF for Monte Carlo (1/Ai)x(1/Aj) */
                const double pdf =
                    (1.0 / patch_area[&recs[i]][ia * recs[i].b_num + ib]) *
                    (1.0 / patch_area[&recs[j]][ja * recs[j].b_num + jb]);

                /* Determine rays of NixNi uniform samples of patch
                   on i to NjxNj uniform samples of patch on j */
                for (int ias = 0; ias < Ni; ias++) {
                  for (int ibs = 0; ibs < Ni; ibs++) {
                    for (int jas = 0; jas < Nj; jas++) {
                      for (int jbs = 0; jbs < Nj; jbs++) {
                        /* Determine sample points xi, xj on both patches */
                        const double u0 = (double)(ias + 0.5) / Ni,
                                     u1 = (double)(ibs + 0.5) / Ni;
                        const double u2 = (double)(jas + 0.5) / Nj,
                                     u3 = (double)(jbs + 0.5) / Nj;

                        const Vector xi =
                            recs[i].p0 +
                            recs[i].edge_a * ((double)(ia + u0) / recs[i].a_num) +
                            recs[i].edge_b * ((double)(ib + u1) / recs[i].b_num);
                        const Vector xj =
                            recs[j].p0 +
                            recs[j].edge_a * ((double)(ja + u2) / recs[j].a_num) +
                            recs[j].edge_b * ((double)(jb + u3) / recs[j].b_num);

                        /* Check for visibility between sample points */
                        const Vector ij = (xj - xi).normalize();

                        double t;
                        int id;
                        Vector normal;
                        if (intersectScene(Ray(xi, ij), &t, &id, &normal) &&
                            id != j) {
                          continue; /* If intersection with other rectangle */
                        }

                        /* Cosines of angles beteen normals and ray inbetween */
                        const double d0 = normal_i.dotProduct(ij);
                        const double d1 = normal_j.dotProduct(-1.0 * ij);

                        /* Continue if patches facing each other */
                        if (d0 > 0.0 && d1 > 0.0) {
                          /* Sample form factor */
                          const double K =
                              d0 * d1 / (M_PI * (xj - xi).lengthSquared());

                          /* Add weighted sample to estimate */
                          F += K / pdf;
                        }
                      }
                    }
                  }
                }

                /* Divide by number of samples */
                F /= (Ni) * (Ni) * (Nj) * (Nj);

                form_factor[&recs[i]][ia * recs[i].b_num + ib][&recs[j]][ja * recs[j].b_num + jb] = F;
                form_factor[&recs[j]][ja * recs[j].b_num + jb][&recs[i]][ia * recs[i].b_num + ib] = F;
              }
            }
          }
        }
      }
    }

    cout << endl;
  }

  /* Divide by area to get final form factors */
  for (auto &rec_a : recs) {
    for (size_t p_a = 0; p_a < rec_a.patch.size(); p_a++) {
      for (auto &rec_b : recs) {
        for (size_t p_b = 0; p_b < rec_b.patch.size(); p_b++) {
          const auto area = patch_area[&rec_a][p_a];
          form_factor[&rec_a][p_a][&rec_b][p_b] = clamp(form_factor[&rec_a][p_a][&rec_b][p_b] / area, 0.0, 1.0);
        }
      }
    }
  }
}

/******************************************************************
 * Iterative computation of radiosity via Gathering; i.e. solution
 * using Gauss-Seidel iteration - reuse already computed values;
 * run-time O(n^2)
 *******************************************************************/

void calculateRadiosity() {
  for (auto &rec_a : recs) {
    for (int ia = 0; ia < rec_a.a_num; ia++) {
      for (int ib = 0; ib < rec_a.b_num; ib++) {
        Color B;

        for (auto &rec_b : recs) {
          for (int ja = 0; ja < rec_b.a_num; ja++) {
            for (int jb = 0; jb < rec_b.b_num; jb++) {
              const double Fij = form_factor[&rec_a][ia * rec_a.b_num + ib][&rec_b][ja * rec_b.b_num + jb];

              /* Add form factor multiplied with radiosity of previous step */
              if (Fij > 0.0)
                B = B + Fij * rec_b.patch[ja * rec_b.b_num + jb];
            }
          }
        }

        /* Multiply sum with color of patch and add emission */
        B = rec_a.color.entrywiseProduct(B) + rec_a.emission;

        /* Store overall patch radiosity of current iteration */
        rec_a.patch[ia * rec_a.b_num + ib] = B;
      }
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

Color radiance(const Ray &ray, bool interpolation = true) {
  double t;
  int id;
  Vector normal;

  /* Find intersected rectangle */
  if (!intersectScene(ray, &t, &id, &normal)) {
    return BackgroundColor;
  }

  /* Determine intersection point on rectangle */
  const Rectangle &obj = recs[id];
  const Vector hitpoint = ray.org + t * ray.dir;

  /* Determine intersected patch */
  const Vector v = hitpoint - obj.p0;
  const double a_len = v.dotProduct(obj.edge_a.normalize());
  const double b_len = v.dotProduct(obj.edge_b.normalize());

  double da = obj.a_num * a_len / obj.a_len;
  double db = obj.b_num * b_len / obj.b_len;

  int ia = int(da);
  if (ia >= obj.a_num)
    ia--;
  int ib = int(db);
  if (ib >= obj.b_num)
    ib--;

  /* Bicubic interpolation for smooth image */
  if (interpolation) {
    Color c[4][4];

    int ia = int(da - 0.5);
    int ib = int(db - 0.5);

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        c[i][j] = obj.sample_patch(ia + i - 1, ib + j - 1);
      }
    }

    int ia0 = int(da - 0.5);
    int ib0 = int(db - 0.5);
    double dx = clamp(da - ia0 - 0.5, 0.0, 1.0);
    double dy = clamp(db - ib0 - 0.5, 0.0, 1.0);

    return ColorUtils::bicubicInterpolate(c, dx, dy) / M_PI;
  } else {
    return obj.patch[ia * obj.b_num + ib] / M_PI;
  }
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

  cout << "Calculating form factors" << endl;
  int patches_a = 12;
  int patches_b = 12;
  int MC_samples = 3;

  calculateFormFactors(patches_a, patches_b, MC_samples);

  /* Iterative solution of radiosity linear system */
  cout << "Calculating radiosity" << endl;
  int iterations = 40;
  for (int i = 0; i < iterations; i++) {
    cout << i << " ";
    calculateRadiosity();
  }
  cout << endl;

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

            /* Determine constant radiance */
            accumulated_radiance =
                accumulated_radiance +
                radiance(Ray(start, dir.normalize()), false) / samples;

            /* Determine interpolated radiance */
            accumulated_radiance2 =
                accumulated_radiance2 +
                radiance(Ray(start, dir.normalize()), true) / samples;
          }

          img.addColor(x, y, accumulated_radiance);
          img_interpolated.addColor(x, y, accumulated_radiance2);
        }
      }
    }
  }

  cout << endl;

  img.save(string("image_patches.ppm"));
  img_interpolated.save(string("image_smooth.ppm"));
}
