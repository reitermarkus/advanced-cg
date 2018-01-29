#ifndef __HSV_H__
#define __HSV_H__

#include <iostream>
#include <cmath>

using namespace std;

struct HSV {
  float h;
  float s;
  float v;

  HSV(const float h_, const float s_, const float v_): h(h_), s(s_), v(v_) {};

  float hue_as_wavelength() const;
  static float randomHue();
  static HSV withRandomHue(float s, float v);
};

ostream& operator<<(std::ostream &strm, const HSV &c);

#endif // __HSV_H__
