#include <iostream>
#include "Sampler.h"

using namespace std;

int main() {
  srand(time(0));

  for (auto i = 0; i < 100; i++) {
    auto randomSpherePoint = Sampler::pointOnSphere();

    cout << randomSpherePoint << endl;
  }

  return 0;
}
