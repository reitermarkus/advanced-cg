#include "Color.h"

#include "macro.h"

#include <random>

const float WAVELENGTH_MAX = 620;
const float WAVELENGTH_RANGE = 170;
const float HUE_RANGE = 270;

float HSV::hueAsWavelength() const {
  return WAVELENGTH_MAX - WAVELENGTH_RANGE / HUE_RANGE * rad_to_deg(this->h);
}

float HSV::randomHue() {
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0.0, 1.0);

  return dis(gen) * HUE_RANGE;
}

HSV HSV::withRandomHue(float s, float v) {
  auto hue = randomHue();
  return HSV(hue, s, v);
}

ostream& operator<<(std::ostream &strm, const HSV &c) {
  return strm << "HSV(" << c.h << ", " << c.s << ", " << c.v << ")";
}

HSV HSV::from(const RGB& color) {
  float max = fmax(fmax(color.r, color.g), color.b);
  float min = fmin(fmin(color.r, color.g), color.b);
  float delta = max - min;

  HSV hsv = HSV(0, 0, max);

  if (delta > 0) {
    if (max > 0) { hsv.s = delta / max; }

    if (max == color.r) { hsv.h = 60.0 * fmod(((color.g - color.b) / delta), 6); } else
    if (max == color.g) { hsv.h = 60.0 * (((color.b - color.r) / delta) + 2.0); } else
    if (max == color.b) { hsv.h = 60.0 * (((color.r - color.g) / delta) + 4.0); }

    while (hsv.h < 0) { hsv.h += 360.0; }
  }

  return hsv;
}
