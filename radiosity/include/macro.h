#ifndef __MACRO_H__
#define __MACRO_H__

#undef M_PI
const double M_PI = atan(1) * 4;

#ifndef drand48
  #define drand48() ((double)rand() / (double)RAND_MAX)
#endif

#ifndef MIN
  #define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#endif

#ifndef MAX
  #define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#endif

#ifndef clamp
  #define clamp(x, min, max) (MIN(max, MAX(x, min)))
#endif

#endif //__MACRO_H__
