#ifndef __MACRO_H__
#define __MACRO_H__

#ifndef drand48
  #define drand48() ((double)rand() / (double)RAND_MAX)
#endif

#undef M_PI
const double M_PI = atan(1) * 4;

double non_uniform_filter_sample() {
  // Get uniform filter sample.
  const double r = 2.0 * drand48();

  // Transform uniform into non-uniform filter sample.
  return r < 1.0 ? (sqrt(r) - 1.0) : (1.0 - sqrt(2.0 - r));
}

#endif // __MACRO_H__
