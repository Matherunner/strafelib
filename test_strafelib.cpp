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

TEST_CASE("friction on velocity", "[friction]") {
    SECTION("geometric friction at 1000fps") {
        double vel[3] = {300, 400, 500};
        fric_vel(vel, 500, 100, 4. / 1000);
        REQUIRE(vel[0] == Approx(298.8));
        REQUIRE(vel[1] == Approx(398.4));
        REQUIRE(vel[2] == 500);
    }
    SECTION("arithmetic friction at 1000fps") {
        double vel[3] = {30, 40, 50};
        fric_vel(vel, 50, 100, 4. / 1000);
        REQUIRE(vel[0] == Approx(29.76));
        REQUIRE(vel[1] == Approx(39.68));
        REQUIRE(vel[2] == Approx(50));
    }
}

TEST_CASE("fme on speed", "[fme]") {
    SECTION("gamma1 at 1000fps") {
        REQUIRE(fme_speed(320, 0.0175, 30, 3.2) == Approx(320.0719919018));
    }
    SECTION("gamma2 at 100fps") {
        REQUIRE(fme_speed(100, 0, 30, 32) == Approx(104.4030650891055));
    }
    SECTION("mu = 0 at 100fps") {
        REQUIRE(fme_speed(1000, 1, 320, 32) == 1000);
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

TEST_CASE("fme minaccel on speed", "[fme]") {
    SECTION("air at 1000fps") {
        REQUIRE(fme_minaccel_speed(2000, 30, 3.2) == 1996.8);
    }
    SECTION("air at 100fps") {
        REQUIRE(fme_minaccel_speed(1200, 30, 32) == 1168);
    }
    SECTION("air at 1000fps, speed = 1") {
        REQUIRE(fme_minaccel_speed(1, 30, 3.2) == 2.2);
    }
    SECTION("air at 100fps, speed = 2") {
        REQUIRE(fme_minaccel_speed(2, 30, 32) == 30);
    }
}

TEST_CASE("fme on velocity", "[fme]") {
    SECTION("120 degrees, air at 1000fps") {
        double vel[2] = {800, 500};
        double speed = std::sqrt(dot_product<2>(vel, vel));
        double theta = 120. * M_PI / 180;
        fme_vel_theta(vel, speed, std::cos(theta), -std::sin(theta), 30, 0.001 * 320 * 10);
        REQUIRE(vel[0] == Approx(797.1744266));
        REQUIRE(vel[1] == Approx(501.5020435));
    }
    SECTION("90 degrees, air at 100fps") {
        double vel[2] = {0, -300};
        double speed = std::sqrt(dot_product<2>(vel, vel));
        fme_vel_theta(vel, speed, 0, 1, 30, 0.01 * 320 * 10);
        REQUIRE(vel[0] == Approx(-30));
        REQUIRE(vel[1] == Approx(-300));
    }
    SECTION("0 degrees, air at 1000fps") {
        double vel[2] = {-500, 500};
        double speed = std::sqrt(dot_product<2>(vel, vel));
        fme_vel_theta(vel, speed, 1, 0, 30, 0.001 * 320 * 10);
        REQUIRE(vel[0] == Approx(-500));
        REQUIRE(vel[1] == Approx(500));
    }
    SECTION("0 degrees, ground at 250fps") {
        double vel[2] = {-100, 0};
        double speed = std::sqrt(dot_product<2>(vel, vel));
        fme_vel_theta(vel, speed, 1, 0, 320, 0.004 * 320 * 10);
        REQUIRE(vel[0] == Approx(-112.8));
        REQUIRE(vel[1] == 0);
    }
}

TEST_CASE("fme_vel_theta benchmark") {
    const double theta = 92 * M_PI / 180;
    const double costheta = std::cos(theta);
    const double sintheta = std::sin(theta);

    BENCHMARK("2000 frames") {
        double v[2] = {80, 50};
        for (int i = 0; i < 2000; ++i) {
            double speed = std::sqrt(dot_product<2>(v, v));
            fme_vel_theta(v, speed, costheta, sintheta, 30, 0.001 * 320 * 100);
        }
        return v[0] + v[1];
    };

    BENCHMARK("2000 frames precomputed speeds theta = 90deg") {
        double v[2] = {80, 50};
        double speedsq = dot_product<2>(v, v);
        double speeds[2000];
        double C = fme_maxaccel_speed_C(speedsq, 30, 32);
        for (int i = 0; i < 2000; ++i) {
            speeds[i] = std::sqrt(speedsq + i * C);
        }
        for (int i = 0; i < 2000; ++i) {
            fme_vel_theta(v, speeds[i], 0, 1, 30, 0.001 * 320 * 100);
        }
        return v[0] + v[1];
    };

    BENCHMARK("2000 frames precomputed speeds theta = zeta") {
        double v[2] = {80, 50};
        double speedsq = dot_product<2>(v, v);
        double speeds[2000];
        double costheta[2000];
        double sintheta[2000];
        double C = fme_maxaccel_speed_C(speedsq, 30, 3.2);
        for (int i = 0; i < 2000; ++i) {
            speeds[i] = std::sqrt(speedsq + i * C);
        }
        for (int i = 0; i < 2000; ++i) {
            fme_maxaccel_cossin_theta(speeds[i], 30, 3.2, costheta + i, sintheta + i);
        }
        for (int i = 0; i < 2000; ++i) {
            fme_vel_theta(v, speeds[i], costheta[i], sintheta[i], 30, 3.2);
        }
        return v[0] + v[1];
    };
}

TEST_CASE("fme maxaccel on speed C", "[fme]") {
    SECTION("air at 1000fps") {
        double speedsq = 1000 * 1000;
        double C = fme_maxaccel_speed_C(speedsq, 30, 3.2);
        double new_speed = std::sqrt(speedsq + C);
        REQUIRE(new_speed == Approx(1000.0908758707881));
    }
    SECTION("air at 1000fps, A = -10") {
        double speedsq = 700 * 700;
        double C = fme_maxaccel_speed_C(speedsq, 30, -3.2);
        double new_speed = std::sqrt(speedsq + C);
        REQUIRE(new_speed == Approx(703.2));
    }
    SECTION("air at 1000fps, zero initial speed") {
        double C = fme_maxaccel_speed_C(0, 30, 3.2);
        double new_speed = std::sqrt(C);
        REQUIRE(new_speed == Approx(3.2));
    }
}

TEST_CASE("collision velocity", "[collision]") {
    SECTION("2D plane") {
        double v[2] = {1000, 0};
        double n[2] = {-3. / 5, 4. / 5};
        collision_vel<2>(v, n, 1);
        REQUIRE(v[0] == Approx(640));
        REQUIRE(v[1] == Approx(480));
    }
}

TEST_CASE("snark_hunt_vel", "[snark]") {
    SECTION("2D vertical") {
        double v[2] = {0, -115};
        double dir[2] = {0, 1};
        snark_hunt_vel<2>(v, dir);
        REQUIRE(v[0] == 0);
        REQUIRE(v[1] == Approx(254));
    }
}

TEST_CASE("tau_g_to_p", "[game]") {
    SECTION("72 fps") {
        REQUIRE(tau_g_to_p(1. / 72) == Approx(0.013));
    }
    SECTION("2000 fps") {
        REQUIRE(tau_g_to_p(1. / 2000) == 0);
    }
    SECTION("1000 fps") {
        REQUIRE(tau_g_to_p(1. / 1000) == Approx(0.001));
    }
    SECTION("100 fps") {
        REQUIRE(tau_g_to_p(1. / 100) == Approx(0.01));
    }
    SECTION("501 fps") {
        REQUIRE(tau_g_to_p(1. / 501) == Approx(0.001));
    }
}

TEST_CASE("water_vel", "[water]") {
    SECTION("1000 fps") {
        double v[3] = {100, 0, 0};
        double a[3] = {1, 0, 0};
        for (int i = 0; i < 1; ++i) {
            double speed = sqrt(dot_product<3>(v, v));
            water_vel(v, speed, a, 1 - 0.001 * 4, 320, 0.001 * 320 * 10);
        }
        REQUIRE(v[0] == 102.16);
        REQUIRE(v[1] == 0);
        REQUIRE(v[2] == 0);
    }
}
