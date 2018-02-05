#include <cmath>

using namespace std;

// https://refractiveindex.info/?shelf=other&book=air&page=Ciddor
inline double refractive_index_of_air(double wavelength_in_nm) {
  double wavelength_in_um = wavelength_in_nm / 1000.0;

  double lambda = wavelength_in_um;
  double lambda_power_minus_2 = pow(lambda, -2.0);

  return 1.0 + 0.05792105 / (238.0185 - lambda_power_minus_2) + 0.00167917 / (57.362 - lambda_power_minus_2);
}
