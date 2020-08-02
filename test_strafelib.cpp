#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "strafelib.hpp"

TEST_CASE("friction on speed", "[friction]") {
    SECTION("geometric friction at 1000fps") {
        REQUIRE(fric_speed(320, 100, 4. / 1000) == Approx(318.72));
    }
    SECTION("arithmetic friction at 1000fps") {
        REQUIRE(fric_speed(80, 100, 4. / 1000) == Approx(79.6));
    }
}

TEST_CASE("fme on speed", "[fme]") {
    SECTION("gamma1 at 1000fps") {
        REQUIRE(fme_speed(320, 0.0175, 30, 3.2) == Approx(320.0719919018));
    }
}

TEST_CASE("fme maxaccel on speed", "[fme]") {
    SECTION("air at 1000fps") {
        REQUIRE(fme_maxaccel_speed(1000, 30, 3.2) == Approx(1000.0908758707881));
    }
    SECTION("air at 100fps") {
        REQUIRE(fme_maxaccel_speed(700, 30, 32) == Approx(700.6425622241344));
    }
    SECTION("air at 1000fps, zero initial speed") {
        REQUIRE(fme_maxaccel_speed(0, 30, 3.2) == Approx(3.2));
    }
    SECTION("air at 1000fps, A = -10") {
        REQUIRE(fme_maxaccel_speed(700, 30, -3.2) == Approx(703.2));
    }
    SECTION("ground at 250fps") {
        REQUIRE(fme_maxaccel_speed(400, 320, 12.8) == Approx(409.91238088157326));
    }
}

TEST_CASE("fme maxaccel on speed C", "[fme]") {
    SECTION("air at 1000fps") {
        double C = fme_maxaccel_speed_C(1000 * 1000, 30, 3.2);
        double new_speed = std::sqrt(1000 * 1000 + C);
        REQUIRE(new_speed == Approx(1000.0908758707881));
    }
    SECTION("air at 1000fps, A = -10") {
        double C = fme_maxaccel_speed_C(700 * 700, 30, -3.2);
        double new_speed = std::sqrt(700 * 700 + C);
        REQUIRE(new_speed == Approx(703.2));
    }
    SECTION("air at 1000fps, zero initial speed") {
        double C = fme_maxaccel_speed_C(0, 30, 3.2);
        double new_speed = std::sqrt(C);
        REQUIRE(new_speed == Approx(3.2));
    }
}
