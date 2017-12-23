#ifndef __MACRO_H__
#define __MACRO_H__

#ifndef drand48
  #define drand48() ((double)rand() / (double)RAND_MAX)
#endif

#undef M_PI
const double M_PI = atan(1) * 4;
#undef M_1_PI
const double M_1_PI = 1 / (atan(1) * 4);

#endif // __MACRO_H__
