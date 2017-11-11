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

SCENARIO("Testing Triangle for Ray intersection.", "[Triangle, Ray]") {
  GIVEN("a right Triangle at (0,0), (0,2), (2,0)") {
    auto triangle = Triangle(Vector(0.0, 0.0, 0.0), Vector(2.0, 0.0, 0.0), Vector(0.0, 2.0, 0.0), Color(), Color());

    WHEN("when a Ray is pointing towards it from the front") {
      auto ray = Ray(Vector(0.1, 0.1, -1.0), Vector(0.0, 0.0, 1.0).normalize());

      THEN("there is an intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection != nullptr);
      }
    }

    WHEN("when a Ray is pointing towards it from behind") {
      auto ray = Ray(Vector(0.1, 0.1, 1.0), Vector(0.0, 0.0, -1.0).normalize());

      THEN("there is an intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection != nullptr);
      }
    }

    WHEN("when a Ray is pointing away from it") {
      auto ray = Ray(Vector(0.1, 0.1, -1.0), Vector(0.0, 0.0, -1.0).normalize());

      THEN("there is no intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection == nullptr);
      }
    }

    WHEN("when a Ray is outside of it") {
      auto ray = Ray(Vector(2.0, 2.0, 1.0), Vector(0.0, 0.0, -1.0).normalize());

      THEN("there is no intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection == nullptr);
      }
    }
  }
}