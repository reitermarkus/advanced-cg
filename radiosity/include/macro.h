#ifndef __MACRO_H__
#define __MACRO_H__

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef drand48
		#define drand48() (rand() / (RAND_MAX + 1.0))
#endif

#ifndef MIN
	#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#endif

#ifndef MAX
	#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#endif

#ifndef clamp
		#define clamp(x, upper, lower) (MIN(upper, MAX(x, lower)))
#endif

#endif //__MACRO_H__
