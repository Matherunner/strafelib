# cppstrafe

Header-only implementations of Half-Life physics primitives with performance as the primary goal. The second goal is the lack of setup and the ease of use. No special tooling or libraries to install before you can start experimenting with Half-Life physics, beyond the standard C++ toolchain you can install with a single command in most macOS and Linux distribution, or with a download of Visual Studio on Windows. The compute routine APIs are all plain and bare.

Just import one header file:

```c++
#include "strafelib.hpp"
```

There's no excuse to not use this library for your research when you need performance!

Most of the compute routines do not stop you from shooting yourself in the foot. For example, the behaviour when passing in negative speeds is undefined.

## Building

When building your research code with this library, it is recommended that you pass in these flags:

    -flto -Ofast -mtune=native -march=native

