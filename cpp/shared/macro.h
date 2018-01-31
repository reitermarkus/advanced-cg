#ifndef __MACRO_H__
#define __MACRO_H__

inline double drand48() {
  return (double)rand() / (double)RAND_MAX;
}

#undef M_PI
const double M_PI = atan(1) * 4.0;

inline double deg_to_rad(double degrees) {
  return degrees * M_PI / 180.0;
}

inline double rad_to_deg(double radians) {
  return radians * 180.0 / M_PI;
}

inline double non_uniform_filter_sample() {
  // Get uniform filter sample.
  const double r = 2.0 * drand48();

  // Transform uniform into non-uniform filter sample.
  return r < 1.0 ? (sqrt(r) - 1.0) : (1.0 - sqrt(2.0 - r));
}

#endif // __MACRO_H__
