#include "ColorUtils.h"

/******************************************************************
 * Helper functions for smooth bicubic (Catmull-Rom) interpolation
 * using 4x4 color patches;
 * First interpolate in y, followed by interpolation of results in x
 *******************************************************************/

Color ColorUtils::bicubicInterpolate(Color p[4][4], double x, double y) {
  Color arr[4];

  auto cubicInterpolate = [] (Color p[4], double x) -> Color {
    return p[1] + 0.5 * x *
                  (p[2] - p[0] +
                    x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] +
                        x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));
  };

  arr[0] = cubicInterpolate(p[0], y);
  arr[1] = cubicInterpolate(p[1], y);
  arr[2] = cubicInterpolate(p[2], y);
  arr[3] = cubicInterpolate(p[3], y);

  return cubicInterpolate(arr, x);
}