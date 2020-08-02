#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "strafelib.hpp"

TEST_CASE("friction on speed", "[friction]") {
    SECTION("geometric friction at 1000 fps, default settings") {
        REQUIRE(fric_speed(320, 100, 4. / 1000) == 318.72);
    }
}
