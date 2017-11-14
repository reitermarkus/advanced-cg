#include "SimpleTriangle.h"

SimpleTriangle::SimpleTriangle(const Vector& a_, const Vector& b_, const Vector& c_) : a(a_), b(b_), c(c_) { }

SimpleTriangle::SimpleTriangle(const Vector &a_, const Vector &b_,
  const Vector &c_, const Color &emission, const Color &color) : Shape(emission, color) {
    a = a_;
    b = a_ + b_;
    c = a_ + c_;
}
