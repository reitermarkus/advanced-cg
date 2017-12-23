#ifndef __MACRO_H__
#define __MACRO_H__

#ifndef drand48
  #define drand48() ((double)rand() / (double)RAND_MAX)
#endif

#endif // __MACRO_H__
