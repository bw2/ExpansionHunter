# Installation

Expansion Hunter is designed for Linux and macOS operating systems. A compiled
binary for the latest release can be downloaded from
[here](https://github.com/Illumina/ExpansionHunter/releases). If you wish to
build the program from source follow the instructions below.

## Building from source

Prerequisites:

 - A recent version of [gcc](https://gcc.gnu.org/) or
   [clang](http://clang.llvm.org/) compiler supporting the C++11 standard
   - The minimum gcc version is 5.1
 - [CMake](https://cmake.org/) version 3.13.0 or above
 - Additional development packages, which depend on the operating system:
     - Centos8
       - `bzip2-devel libcurl-devel libstdc++-static openssl-devel xz-devel zlib-devel`
     - Ubuntu 20.04
       - `libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev zlib1g-dev `
     - macOS 10.15
       - `xz` (from homebrew)

If the above prerequisites are satisfied, you are ready to
build the program. Note that during the build procedure, cmake will
attempt to download and install `abseil`, `boost`, `googletest`, `htslib`,
and `spdlog` so an active internet connection is required. Assuming
that the source code is contained in a directory `ExpansionHunter/`,
the build procedure can be initiated as follows:

```bash
$ cd ExpansionHunter
$ mkdir build
$ cd build
$ cmake ..
$ make
```

If all the above steps were successful, the ExpansionHunter executable can be found in:

    build/install/bin/ExpansionHunter

### Building without internet access

Each dependency's source archive can be given as a local file instead of a URL. Download the same files on a
machine with internet access (the default URLs are listed in the top-level `CMakeLists.txt`), copy them over, and
pass their paths:

```bash
$ cmake .. \
    -DLIBDEFLATE_SOURCE_ARCHIVE=/path/to/libdeflate-1.23.tar.gz \
    -DHTSLIB_SOURCE_ARCHIVE=/path/to/htslib-1.21.tar.bz2 \
    -DBOOST_SOURCE_ARCHIVE=/path/to/boost_1_87_0.tar.gz \
    -DSPDLOG_SOURCE_ARCHIVE=/path/to/spdlog-1.15.1.tar.gz \
    -DGOOGLETEST_SOURCE_ARCHIVE=/path/to/googletest-1.16.0.tar.gz
```

### Building against installed htslib or Boost

By default the build downloads and compiles its own htslib (with libdeflate) and Boost. To use copies already
installed on your machine instead:

```bash
$ cmake .. \
    -DUSE_SYSTEM_HTSLIB=ON -DSYSTEM_HTSLIB_PREFIX=/path/to/htslib \
    -DUSE_SYSTEM_BOOST=ON -DSYSTEM_BOOST_PREFIX=/path/to/boost
```

 - `SYSTEM_HTSLIB_PREFIX` and `SYSTEM_BOOST_PREFIX` are optional. Without them, the standard install locations
   are searched. Use a fresh build directory when switching these options.
 - htslib must be version 1.10 or newer and is linked as a shared library. Add
   `-DLINK_SYSTEM_HTSLIB_STATICALLY=ON` to link its static `libhts.a` instead; libdeflate is then linked too
   (a static `libdeflate.a` is preferred, searched for in the htslib prefix first). Reading `gs://`, `s3://` or
   `https://` input requires an htslib built with libcurl support (`--enable-libcurl`, plus `--enable-gcs` and
   `--enable-s3` for those schemes).
 - Boost must be version 1.84 or newer and include static libraries for `program_options`, `filesystem`,
   `system` and `iostreams`.
 - The installed ExpansionHunter binary does not record where a shared htslib is. If it lives outside the standard
   library locations, add its `lib` directory to `LD_LIBRARY_PATH` when running ExpansionHunter on Linux. On
   macOS the full path of the htslib library is recorded in the binary by the system linker.

A binary built this way depends on the exact libraries installed on the machine that built it. Use it on that
machine only: do not distribute it, and when reporting a problem, mention which htslib version it was built
against. The Docker image is always built with the default settings.
