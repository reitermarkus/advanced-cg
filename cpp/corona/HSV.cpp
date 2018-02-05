#include "HSV.h"

#include "macro.h"
#include "Vector.h"

#include <random>

const float WAVELENGTH_MAX = 620;
const float WAVELENGTH_RANGE = 170;
const float HUE_RANGE = 270;

float HSV::hueAsWavelength() const {
  return WAVELENGTH_MAX - WAVELENGTH_RANGE / HUE_RANGE * this->h;
}

float HSV::randomHue() {
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0.0, HUE_RANGE);

  return dis(gen);
}

HSV HSV::withRandomHue(float s, float v) {
  auto hue = randomHue();
  return HSV(hue, s, v);
}

ostream& operator<<(std::ostream &strm, const HSV &c) {
  return strm << "HSV(" << c.h << ", " << c.s << ", " << c.v << ")";
}

HSV HSV::rgbToHsv(const Color& color) {
  float max = fmax(fmax(color.x, color.y), color.z);
  float min = fmin(fmin(color.x, color.y), color.z);
  float delta = max - min;

  HSV hsv = HSV(0, 0, max / 255);

  if (delta > 0) {
    if (max > 0) { hsv.s = delta / max; }

    if (max == color.x) { hsv.h = 60.0 * fmod(((color.y - color.z) / delta), 6); } else
    if (max == color.y) { hsv.h = 60.0 * (((color.z - color.x) / delta) + 2.0); } else
    if (max == color.z) { hsv.h = 60.0 * (((color.x - color.y) / delta) + 4.0); }

    while (hsv.h < 0) { hsv.h += 360.0; }
  }

  return hsv;
}

Color HSV::toRGB() const {
  float hue = this->h;

  while (hue < 0) { hue += 360; }

  float chroma = this->v * this->s;
  float hue_mod = fmod(hue / 60.0, 6);
  float m = this->v - chroma;
  float x = chroma * (1.0 - fabs(1.0 - fmod(hue_mod, 2))) + m;

  switch ((int)hue_mod) {
    case 0: return Color(this->v, x,       m);
    case 1: return Color(x,       this->v, m);
    case 2: return Color(m,       this->v, x);
    case 3: return Color(m,       x,       this->v);
    case 4: return Color(x,       m,       this->v);
    case 5: return Color(this->v, m,       x);
  }
}
