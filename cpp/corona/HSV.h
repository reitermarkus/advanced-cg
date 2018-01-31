#ifndef __HSV_H__
#define __HSV_H__

#include <iostream>
#include <cmath>

#include "../shared/Vector.h"

using namespace std;

struct HSV {
  float h;
  float s;
  float v;

  HSV(const float h_, const float s_, const float v_): h(h_), s(s_), v(v_) {};

  float hueAsWavelength() const;
  static float randomHue();
  static HSV withRandomHue(float s, float v);
  static HSV rgbToHsv(const Color& color);
};

ostream& operator<<(std::ostream &strm, const HSV &c);

#endif // __HSV_H__
