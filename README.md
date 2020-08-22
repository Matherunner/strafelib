# strafelib (WIP!)

[![Build Status](https://travis-ci.com/Matherunner/strafelib.svg?branch=master)](https://travis-ci.com/Matherunner/strafelib)

**NOTE**: This library is work-in-progress!

Header-only implementations of Half-Life physics primitives with performance as the primary goal. The secondary goal is the minimal setup and easy of use. No special tooling or libraries to hunt down and install before you can start experimenting with Half-Life physics, beyond the standard C++ toolchain you can install with a single command in most macOS and Linux distribution, or with a download of Visual Studio on Windows. The compute routine APIs are all plain and bare.

## How to use

Just import one header file:

```c++
#include "strafelib.hpp"
```

When building your research code with this library, it is recommended that you pass in these flags to GCC:

    -flto -Ofast -mtune=native -march=native

## Performance

I will give you an idea of the single-core performance of this library. My CPU is a stock [Intel Core i7-8700](https://ark.intel.com/content/www/us/en/ark/products/126686/intel-core-i7-8700-processor-12m-cache-up-to-4-60-ghz.html).

- Using the `fme_vel_theta` with default Half-Life settings and 1000 fps, without any precomputations, it calculates roughly 80,000,000 velocity vectors per second.
- By precomputing speeds using `fme_maxaccel_speed_C` at 100aa and 1000 fps, and using the speeds in `fme_vel_theta`, it calculates roughly 500,000,000 velocity vectors per second.

## Building

No building is needed to start using this library in your own code. However, you may wish to build the tests for this library by running

    $ make test
    $ ./test_strafelib
