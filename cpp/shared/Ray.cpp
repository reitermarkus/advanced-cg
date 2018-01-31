#include "Ray.h"

#include <iostream>

ostream& operator<<(std::ostream &strm, const Ray &r) {
  return strm << "Ray(" << r.org << ", " << r.dir << ")";
}
