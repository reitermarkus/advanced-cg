#ifndef __DIFFRACTION_H__
#define __DIFFRACTION_H__

#include <cmath>

// https://en.wikipedia.org/wiki/Diffraction
//   1. Single-slit diffraction
//   2. Diffraction-limited imaging
inline double minDiffractionAngle(double wavelength, double aperture, double focal_length) {
  auto f_number = focal_length / aperture;
  auto d = 1.22 * wavelength * f_number; // 2

  auto sin_angle = wavelength / d; // 1

  return asin(sin_angle);
}

inline double diffractionAngle(double wavelength) {
  return wavelength / 270.0 * deg_to_rad(45);
}

#endif // __DIFFRACTION_H__
