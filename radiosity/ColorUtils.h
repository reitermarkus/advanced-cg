#ifndef __COLOR_UTILS_H__
#define __COLOR_UTILS_H__

#include "Vector.h"

struct ColorUtils {
    static Color bicubicInterpolate(Color p[4][4], double x, double y);
};

#endif //__COLOR_UTILS_H__