#include "HSV.h"

#include "macro.h"

#include <random>

const float WAVELENGTH_MAX = 620;
const float WAVELENGTH_RANGE = 170;
const float HUE_RANGE = 270;

float HSV::hue_as_wavelength() const {
  return WAVELENGTH_MAX - WAVELENGTH_RANGE / HUE_RANGE * rad_to_deg(this->h);
}

float HSV::random_hue() {
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0.0, 1.0);

  return dis(gen) * HUE_RANGE;
}

HSV HSV::with_random_hue(float s, float v) {
  auto hue = random_hue();
  return HSV(hue, s, v);
}

ostream& operator<<(std::ostream &strm, const HSV &c) {
  return strm << "HSV(" << c.h << ", " << c.s << ", " << c.v << ")";
}
