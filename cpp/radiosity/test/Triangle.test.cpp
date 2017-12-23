#include "catch.hpp"

#include "Triangle.h"

SCENARIO("Calculating Triangle patch count.", "[Triangle]") {
  GIVEN("any Triangle") {
    auto triangle = Triangle(Vector(), Vector(), Vector(), Color(), Color());

    WHEN("the sides are divided into 1 part") {
      triangle.init_patches(1);

      THEN("the patch count is 1") {
        REQUIRE(triangle.patches.size() == 1);
      }
    }

    WHEN("the sides are divided into 2 parts") {
      triangle.init_patches(2);

      THEN("the patch count is 4") {
        REQUIRE(triangle.patches.size() == 4);
      }
    }


    WHEN("the sides are divided into 3 parts") {
      triangle.init_patches(3);

      THEN("the patch count is 9") {
        REQUIRE(triangle.patches.size() == 9);
      }
    }

    WHEN("the sides are divided into 4 parts") {
      triangle.init_patches(4);

      THEN("the patch count is 16") {
        REQUIRE(triangle.patches.size() == 16);
      }
    }
  }
}

SCENARIO("Testing Triangle for Ray intersection.", "[Triangle, Ray]") {
  GIVEN("a right Triangle at (2,2), (-2,0), (0,-2)") {
    auto triangle = Triangle(Vector(2.0, 2.0, 0.0), Vector(2.0, 0.0, 0.0), Vector(0.0, 2.0, 0.0), Color(), Color());

    WHEN("when a Ray is pointing towards it from the front") {
      auto ray = Ray(Vector(2.25, 2.25, -1.0), Vector(0.0, 0.0, 1.0).normalize());

      THEN("there is an intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection > 0.0);
      }
    }
  }

  GIVEN("a right Triangle at (0,0), (2,0), (0,2)") {
    auto triangle = Triangle(Vector(0.0, 0.0, 0.0), Vector(2.0, 0.0, 0.0), Vector(0.0, 2.0, 0.0), Color(), Color());

    WHEN("when a Ray is pointing towards it from the front") {
      auto ray = Ray(Vector(0.1, 0.1, -1.0), Vector(0.0, 0.0, 1.0).normalize());

      THEN("there is an intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection > 0.0);
      }
    }

    WHEN("when a Ray is pointing towards it from behind") {
      auto ray = Ray(Vector(0.1, 0.1, 1.0), Vector(0.0, 0.0, -1.0).normalize());

      THEN("there is an intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection > 0.0);
      }
    }

    WHEN("when a Ray is pointing away from it") {
      auto ray = Ray(Vector(0.1, 0.1, -1.0), Vector(0.0, 0.0, -1.0).normalize());

      THEN("there is no intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection == 0.0);
      }
    }

    WHEN("when a Ray is outside of it") {
      auto ray = Ray(Vector(2.0, 2.0, 1.0), Vector(0.0, 0.0, -1.0).normalize());

      THEN("there is no intersection") {
        auto intersection = triangle.intersect(ray);
        REQUIRE(intersection == 0.0);
      }
    }
  }
}


SCENARIO("Testing Triangle division.", "[Triangle]") {
  GIVEN("a right Triangle at (0,0), (1,0), (0,1)") {
    auto triangle = Triangle(Vector(0.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0), Color(), Color());

    WHEN("divided by 3") {
      THEN("it should output all coordinates") {
        triangle.divide(3);
      }
    }
  }
}


SCENARIO("Testing Triangle sample point.", "[Triangle]") {
  GIVEN("a right Triangle at (0,0), (1,0), (0,1)") {
    auto p0 = Vector(0.0, 0.0, 0.0);
    auto p1 = Vector(1.0, 0.0, 0.0);
    auto p2 = Vector(0.0, 1.0, 0.0);

    auto triangle = Triangle(p0, p1, p2, Color(), Color());

    WHEN("generating a random sample point") {
      THEN("it is inside the triangle") {
        auto q = Triangle::sample(p0, p1, p2);
        auto intersection = triangle.intersect(Ray(q + Vector(0.0, 0.0, -1), Vector(0.0, 0.0, 1.0)));
        REQUIRE(intersection > 0.0);
      }
    }
  }
}
