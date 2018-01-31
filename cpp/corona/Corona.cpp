#include <iostream>
#include "Sampler.h"
#include "Color.h"

using namespace std;

int main() {
  srand(time(0));

  for (auto i = 0; i < 100; i++) {
    auto randomSpherePoint = Sampler::pointOnSphere();

    cout << randomSpherePoint << endl;
  }

  auto hsv = HSV::rgbToHsv(RGB(255, 138, 0));

  cout << hsv.h << " " << hsv.s << " " << hsv.v << endl;


  cout << Sampler::sphericalRay(Vector(0.0, 0.0, 0.0), 2.5) << endl;

  return 0;
}
