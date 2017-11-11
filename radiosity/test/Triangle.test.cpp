#include "catch.hpp"

#include "Triangle.h"

SCENARIO("Calculating Triangle patch count.", "[Triangle]") {
  GIVEN("any Triangle") {
    auto triangle = Triangle(Vector(), Vector(), Vector(), Color(), Color());

    WHEN("the sides are divided into 1 part") {
      triangle.init_patches(1);

      THEN("the patch count is 1") {
        REQUIRE(triangle.patch.size() == 1);
      }
    }

    WHEN("the sides are divided into 2 parts") {
      triangle.init_patches(2);

      THEN("the patch count is 4") {
        REQUIRE(triangle.patch.size() == 4);
      }
    }


    WHEN("the sides are divided into 3 parts") {
      triangle.init_patches(3);

      THEN("the patch count is 9") {
        REQUIRE(triangle.patch.size() == 9);
      }
    }

    WHEN("the sides are divided into 4 parts") {
      triangle.init_patches(4);

      THEN("the patch count is 16") {
        REQUIRE(triangle.patch.size() == 16);
      }
    }
  }
}
