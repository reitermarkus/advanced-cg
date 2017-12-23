#include "../shared/test/catch.hpp"

#include "Vector.h"

SCENARIO("Testing Vector equality.", "[Vector]") {
  GIVEN("two Vector with x: 8, y: 4, z: 2") {
    auto vector_a = Vector(8, 4, 2);
    auto vector_b = Vector(8, 4, 2);
    WHEN("compared with each other") {
      THEN("it should return true") {
        REQUIRE(vector_a == vector_b);
      }
    }
  }

  GIVEN("two Vector with x: 16.65, y: 4.75, z: 9.001") {
    auto vector_a = Vector(16.65, 4.75, 9.001);
    auto vector_b = Vector(16.65, 4.75, 9.001);
    WHEN("compared with each other") {
      THEN("it should return true") {
        REQUIRE(vector_a == vector_b);
      }
    }
  }

  GIVEN("a Vector with x: 9, y: 3.75, z: 2 and a Vector with x: 9, y: 3.7502, z: 2") {
    auto vector_a = Vector(9, 3.75, 2);
    auto vector_b = Vector(9, 3.7502, 2);
    WHEN("compared with each other") {
      THEN("it should return false") {
        REQUIRE_FALSE(vector_a == vector_b);
      }
    }
  }
}
