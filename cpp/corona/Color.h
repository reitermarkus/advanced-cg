#ifndef __COLOR_H__
#define __COLOR_H__

#include <iostream>
#include <cmath>

using namespace std;

struct RGB {
  float r;
  float g;
  float b;

  RGB(const float r_, const float g_, const float b_): r(r_), g(g_), b(b_) {};
};

struct HSV {
  float h;
  float s;
  float v;

  HSV(const float h_, const float s_, const float v_): h(h_), s(s_), v(v_) {};

  float hueAsWavelength() const;
  static float randomHue();
  static HSV withRandomHue(float s, float v);
  static HSV rgbToHsv(const RGB& color);
};

ostream& operator<<(std::ostream &strm, const HSV &c);

#endif // __COLOR_H__
