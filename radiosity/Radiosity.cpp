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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "Vector.h"
#include "Image.h"
#include "Rectangle.h"
#include "Ray.h"

using namespace std;

const double Over_M_PI = 1.0 / M_PI;

static double *form_factor;
static int patch_num = 0;

const Color BackgroundColor(0.0, 0.0, 0.0);

/******************************************************************
 * Hard-coded scene definition: the geometry is composed of rectangles.
 * These are defined by:
 * - vector to corner(origin), edge a, edge b
 * - emitted light energy (light sources), surface reflectivity (~color)
 *******************************************************************/
vector<Rectangle> recs = {
  /* Cornell Box walls */
  Rectangle(Vector(0.0, 0.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 80.0, 0.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Back */
  Rectangle(Vector(0.0, 0.0, 170.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, -170.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Bottom */
  Rectangle(Vector(0.0, 80.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Top */
  Rectangle(Vector(0.0, 0.0, 170.0), Vector(0.0, 0.0, -170.0), Vector(0.0, 80.0, 0.0),
            Color(), Color(0.75, 0.25, 0.25)), /* Left */
  Rectangle(Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0), Vector(0.0, 80.0, 0.0),
            Color(), Color(0.25, 0.25, 0.75)), /* Right */
  Rectangle(Vector(100.0, 0.0, 170.0), Vector(-100.0, 0.0, 0.0), Vector(0.0, -80.0, 0.0),
            Color(), Color(0, 1, 0)), /* Front (not visible) */

  /* Area light source on top */
  Rectangle(Vector(40.0, 79.99, 65.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, 20.0),
            Color(12, 12, 12), Color(0.75, 0.75, 0.75)),

  /* Cuboid in room */
  Rectangle(Vector(30.0, 0.0, 100.0), Vector(0.0, 0.0, -20.0), Vector(0.0, 40.0, 0.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Right */
  Rectangle(Vector(10.0, 0.0, 80.0), Vector(0.0, 0.0, 20.0), Vector(0.0, 40.0, 0.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Left */
  Rectangle(Vector(10.0, 0.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 40.0, 0.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Front */
  Rectangle(Vector(30.0, 0.0, 80.0), Vector(-20.0, 0.0, 0.0), Vector(0.0, -40.0, 0.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Back */
  Rectangle(Vector(10.0, 40.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, -20.0),
            Color(), Color(0.75, 0.75, 0.75)), /* Top */
};

/******************************************************************
 * Check for closest intersection of a ray with the scene;
 * Returns true if intersection is found, as well as ray parameter
 * of intersection and id of intersected object
 *******************************************************************/
bool Intersect_Scene(const Ray &ray, double *t, int *id, Vector *normal) {
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
void Calculate_Form_Factors(const int a_div_num, const int b_div_num,
                            const int mc_sample) {
  /* Total number of patches in scene */
  const int n = recs.size();
  for (int i = 0; i < n; i++) {
    recs[i].init_patchs(a_div_num, b_div_num);
    patch_num += recs[i].a_num * recs[i].b_num;
  }

  cout << "Number of rectangles: " << n << endl;
  cout << "Number of patches: " << patch_num << endl;
  int form_factor_num = patch_num * patch_num;
  cout << "Number of form factors: " << form_factor_num << endl;

  /* 1D-array to hold form factor pairs */
  form_factor = new double[form_factor_num];
  memset(form_factor, 0.0, sizeof(double) * form_factor_num);

  /* 1D-array with patch areas */
  double *patch_area = new double[patch_num];
  memset(patch_area, 0.0, sizeof(double) * patch_num);

  /* Precompute patch areas, assuming same size for each rectangle */
  for (int i = 0; i < n; i++) {
    int patch_i = 0;
    for (int k = 0; k < i; k++)
      patch_i += recs[k].a_num * recs[k].b_num;

    for (int ia = 0; ia < recs[i].a_num; ia++) {
      for (int ib = 0; ib < recs[i].b_num; ib++) {
        patch_area[patch_i + ia * recs[i].b_num + ib] =
            ((recs[i].edge_a / recs[i].a_num)
                 .Cross((recs[i].edge_b / recs[i].b_num)))
                .Length();
      }
    }
  }

  /* Offsets for indexing of patches in 1D-array */
  int *offset = new int[n];

  for (int i = 0; i < n; i++) {
    offset[i] = 0;
    for (int k = 0; k < i; k++)
      offset[i] += recs[k].a_num * recs[k].b_num;
  }

  /* Loop over all rectangles in scene */
  for (int i = 0; i < n; i++) {
    int patch_i = offset[i];

    cout << i << " ";

    /* Loop over all patches in rectangle i */
    for (int ia = 0; ia < recs[i].a_num; ia++) {
      cout << "*" << flush;
      for (int ib = 0; ib < recs[i].b_num; ib++) {
        const Vector normal_i = recs[i].normal;

        int patch_j = 0;

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
                    (1.0 / patch_area[offset[i] + ia * recs[i].b_num + ib]) *
                    (1.0 / patch_area[offset[j] + ja * recs[j].b_num + jb]);

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
                            recs[i].edge_a *
                                ((double)(ia + u0) / recs[i].a_num) +
                            recs[i].edge_b *
                                ((double)(ib + u1) / recs[i].b_num);
                        const Vector xj =
                            recs[j].p0 +
                            recs[j].edge_a *
                                ((double)(ja + u2) / recs[j].a_num) +
                            recs[j].edge_b *
                                ((double)(jb + u3) / recs[j].b_num);

                        /* Check for visibility between sample points */
                        const Vector ij = (xj - xi).Normalized();

                        double t;
                        int id;
                        Vector normal;
                        if (Intersect_Scene(Ray(xi, ij), &t, &id, &normal) &&
                            id != j) {
                          continue; /* If intersection with other rectangle */
                        }

                        /* Cosines of angles beteen normals and ray inbetween */
                        const double d0 = normal_i.Dot(ij);
                        const double d1 = normal_j.Dot(-1.0 * ij);

                        /* Continue if patches facing each other */
                        if (d0 > 0.0 && d1 > 0.0) {
                          /* Sample form factor */
                          const double K =
                              d0 * d1 / (M_PI * (xj - xi).LengthSquared());

                          /* Add weighted sample to estimate */
                          F += K / pdf;
                        }
                      }
                    }
                  }
                }

                /* Divide by number of samples */
                F /= (Ni) * (Ni) * (Nj) * (Nj);

                form_factor[patch_i * patch_num + patch_j] = F;
              }
              patch_j++;
            }
          }
        }
        patch_i++;
      }
    }

    cout << endl;
  }

  /* Copy upper to lower triangular values */
  for (int i = 0; i < patch_num - 1; i++) {
    for (int j = i + 1; j < patch_num; j++) {
      form_factor[j * patch_num + i] = form_factor[i * patch_num + j];
    }
  }

  /* Divide by area to get final form factors */
  for (int i = 0; i < patch_num; i++) {
    for (int j = 0; j < patch_num; j++) {
      form_factor[i * patch_num + j] /= patch_area[i];

      /* Clamp to [0,1] */
      if (form_factor[i * patch_num + j] > 1.0)
        form_factor[i * patch_num + j] = 1.0;
    }
  }
}

/******************************************************************
 * Iterative computation of radiosity via Gathering; i.e. solution
 * using Gauss-Seidel iteration - reuse already computed values;
 * run-time O(n^2)
 *******************************************************************/

void Calculate_Radiosity(const int iteration) {
  const int n = recs.size();
  int patch_i = 0;

  for (int i = 0; i < n; i++) {
    for (int ia = 0; ia < recs[i].a_num; ia++) {
      for (int ib = 0; ib < recs[i].b_num; ib++) {
        Color B;

        int patch_j = 0;
        for (int j = 0; j < n; j++) {
          for (int ja = 0; ja < recs[j].a_num; ja++) {
            for (int jb = 0; jb < recs[j].b_num; jb++) {
              const double Fij = form_factor[patch_i * patch_num + patch_j];

              /* Add form factor multiplied with radiosity of previous step */
              if (Fij > 0.0)
                B = B + Fij * recs[j].patch[ja * recs[j].b_num + jb];

              patch_j++;
            }
          }
        }
        /* Multiply sum with color of patch and add emission */
        B = recs[i].color.MultComponents(B) + recs[i].emission;

        /* Store overall patch radiosity of current iteration */
        recs[i].patch[ia * recs[i].b_num + ib] = B;
        patch_i++;
      }
    }
  }
}

/******************************************************************
 * Helper functions for smooth bicubic (Catmull-Rom) interpolation
 * using 4x4 color patches;
 * First interpolate in y, followed by interpolation of results in x
 *******************************************************************/

Color cubicInterpolate(Color p[4], double x) {
  return p[1] + 0.5 * x *
                    (p[2] - p[0] +
                     x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
                          x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
}

Color bicubicInterpolate(Color p[4][4], double x, double y) {
  Color arr[4];

  arr[0] = cubicInterpolate(p[0], y);
  arr[1] = cubicInterpolate(p[1], y);
  arr[2] = cubicInterpolate(p[2], y);
  arr[3] = cubicInterpolate(p[3], y);

  return cubicInterpolate(arr, x);
}

/******************************************************************
 * Compute radiance from radiosity by shooting rays into the scene;
 * Radiance directly proportional to radiosity for assumed diffuse
 * emitters/surfaces (multiply by PI);
 * At intersections either constant patch color is returned or a
 * smoothly interpolated color of 4x4 neighboring patches
 *******************************************************************/

Color Radiance(const Ray &ray, const int depth, bool interpolation = true) {
  double t;
  int id;
  Vector normal;

  /* Find intersected rectangle */
  if (!Intersect_Scene(ray, &t, &id, &normal)) {
    return BackgroundColor;
  }

  /* Determine intersection point on rectangle */
  const Rectangle &obj = recs[id];
  const Vector hitpoint = ray.org + t * ray.dir;

  /* Determine intersected patch */
  const Vector v = hitpoint - obj.p0;
  const double a_len = v.Dot(obj.edge_a.Normalized());
  const double b_len = v.Dot(obj.edge_b.Normalized());

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
    double dx = da - ia0 - 0.5;
    double dy = db - ib0 - 0.5;

    if (dx < 0.0)
      dx = 0.0;
    if (dx >= 1.0)
      dx = 1.0;
    if (dy < 0.0)
      dy = 0.0;
    if (dy >= 1.0)
      dy = 1.0;

    return bicubicInterpolate(c, dx, dy) * Over_M_PI;
  } else {
    return obj.patch[ia * obj.b_num + ib] * Over_M_PI;
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

int main(int argc, char **argv) {
  int width = 640;
  int height = 480;
  int samples = 4;

  /* Set camera origin and viewing direction (negative z direction) */
  Ray camera(Vector(50.0, 52.0, 295.6),
             Vector(0.0, -0.042612, -1.0).Normalized());

  /* Image edge vectors for pixel sampling */
  Vector cx = Vector(width * 0.5135 / height);
  Vector cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;

  /* Two final renderings; one with constant, one with interpolated patch colors
   */
  Image img(width, height);
  Image img_interpolated(width, height);

  cout << "Calculating form factors" << endl;
  int patches_a = 12;
  int patches_b = 12;
  int MC_samples = 3;

  Calculate_Form_Factors(patches_a, patches_b, MC_samples);

  /* Iterative solution of radiosity linear system */
  cout << "Calculating radiosity" << endl;
  int iterations = 40;
  for (int i = 0; i < iterations; i++) {
    cout << i << " ";
    Calculate_Radiosity(i);
  }
  cout << endl;

  /* Loop over image rows */
  for (int y = 0; y < height; y++) {
    cout << "\rRendering (" << samples * 4 << " spp) "
         << (100.0 * y / (height - 1)) << "%     ";
    srand(y * y * y);

    /* Loop over row pixels */
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
            const double r1 = 2.0 * drand48();
            const double r2 = 2.0 * drand48();

            /* Transform uniform into non-uniform filter samples */
            double dx = r1 < 1.0 ? (sqrt(r1) - 1.0) : (1.0 - sqrt(2.0 - r1));
            double dy = r2 < 1.0 ? (sqrt(r2) - 1.0) : (1.0 - sqrt(2.0 - r2));

            /* Ray direction into scene from camera through sample */
            Vector dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
                         cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) +
                         camera.dir;

            /* Extend camera ray to start inside box */
            Vector start = camera.org + dir * 130.0;

            /* Determine constant radiance */
            accumulated_radiance =
                accumulated_radiance +
                Radiance(Ray(start, dir.Normalized()), 0, false) / samples;

            /* Determine interpolated radiance */
            accumulated_radiance2 =
                accumulated_radiance2 +
                Radiance(Ray(start, dir.Normalized()), 0, true) / samples;
          }

          img.addColor(x, y, accumulated_radiance);
          img_interpolated.addColor(x, y, accumulated_radiance2);
        }
      }
    }
  }

  cout << endl;

  img.Save(string("image_patches.ppm"));
  img_interpolated.Save(string("image_smooth.ppm"));
}
