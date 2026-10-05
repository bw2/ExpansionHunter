cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check for required packages:
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)

include(GNUInstallDirs)

if(NOT APPLE)
  find_package(OpenSSL REQUIRED)
endif()

include(ExternalProject)
set(installDir ${CMAKE_CURRENT_BINARY_DIR}/install)

# By default htslib (with libdeflate) and Boost are downloaded and built here. These options build against
# copies already installed on this machine instead. Such builds are for use on this machine only; see
# docs/02_Installation.md.
option(USE_SYSTEM_HTSLIB "Use an installed htslib (1.10 or newer) instead of downloading htslib and libdeflate" OFF)
set(SYSTEM_HTSLIB_PREFIX "" CACHE PATH "Install prefix of the htslib used with USE_SYSTEM_HTSLIB (default: search the standard locations)")
option(LINK_SYSTEM_HTSLIB_STATICALLY "With USE_SYSTEM_HTSLIB, link the installed static libhts.a (plus libdeflate) instead of the shared library" OFF)
option(USE_SYSTEM_BOOST "Use an installed Boost (1.84 or newer, static libraries) instead of downloading Boost" OFF)
set(SYSTEM_BOOST_PREFIX "" CACHE PATH "Install prefix of the Boost used with USE_SYSTEM_BOOST (default: search the standard locations)")

# Where to get each dependency's source archive: a URL or a path to an already-downloaded copy of the same file,
# for building without internet access.
set(LIBDEFLATE_SOURCE_ARCHIVE https://github.com/ebiggers/libdeflate/releases/download/v1.23/libdeflate-1.23.tar.gz
	CACHE STRING "URL or local path of the libdeflate source archive")
set(HTSLIB_SOURCE_ARCHIVE https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
	CACHE STRING "URL or local path of the htslib source archive")
set(BOOST_SOURCE_ARCHIVE https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz
	CACHE STRING "URL or local path of the Boost source archive")
set(SPDLOG_SOURCE_ARCHIVE https://github.com/gabime/spdlog/archive/refs/tags/v1.15.1.tar.gz
	CACHE STRING "URL or local path of the spdlog source archive")
set(GOOGLETEST_SOURCE_ARCHIVE https://github.com/google/googletest/archive/refs/tags/v1.16.0.tar.gz
	CACHE STRING "URL or local path of the googletest source archive")

set(ehunterDependencies spdlog googletest)

if(NOT USE_SYSTEM_HTSLIB)
list(APPEND ehunterDependencies htslib)

ExternalProject_Add(libdeflate
	URL ${LIBDEFLATE_SOURCE_ARCHIVE}
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
	-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	-DCMAKE_INSTALL_PREFIX=${installDir}
	-DCMAKE_C_FLAGS=-O2
	-DLIBDEFLATE_BUILD_SHARED_LIB=OFF
	-DLIBDEFLATE_BUILD_STATIC_LIB=ON
	-DLIBDEFLATE_BUILD_GZIP=OFF
	-DCMAKE_POSITION_INDEPENDENT_CODE=ON
)

ExternalProject_Add(htslib
	BUILD_IN_SOURCE YES
	DEPENDS libdeflate
	URL ${HTSLIB_SOURCE_ARCHIVE}
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CONFIGURE_COMMAND env CFLAGS=-O2\ -g\ -I${installDir}/include LDFLAGS=-L${installDir}/${CMAKE_INSTALL_LIBDIR} ./configure --prefix=${installDir} --enable-libcurl --enable-s3 --enable-gcs --with-libdeflate
	BUILD_COMMAND $(MAKE) lib-static
)
endif()


if(NOT USE_SYSTEM_BOOST)
list(APPEND ehunterDependencies Boost)

# Setup boost user config so that it uses the same c++ compiler as other components:
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(BOOST_COMPILER_TAG "gcc")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	set(BOOST_COMPILER_TAG "clang")
	set(B2_OPTIONS ${B2_OPTIONS} "toolset=clang")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	set(BOOST_COMPILER_TAG "darwin")
	set(B2_OPTIONS ${B2_OPTIONS} "toolset=clang")
endif()

if (BOOST_COMPILER_TAG)
	set(BOOST_UCONFIG "${CMAKE_BINARY_DIR}/user-config.jam")
	file(WRITE "${BOOST_UCONFIG}" "using ${BOOST_COMPILER_TAG} : : \"${CMAKE_CXX_COMPILER}\" ;\n")
	set(BOOST_PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${BOOST_UCONFIG} tools/build/src/user-config.jam)
endif()

ExternalProject_Add(Boost
	BUILD_IN_SOURCE YES
	URL ${BOOST_SOURCE_ARCHIVE}
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	PATCH_COMMAND ${BOOST_PATCH_COMMAND}
	CONFIGURE_COMMAND ./bootstrap.sh --prefix=${installDir}
	BUILD_COMMAND ./b2 install -j8 cxxstd=17 variant=release optimization=speed -sNO_ZSTD=1 --with-filesystem --with-system --with-iostreams --with-program_options link=static ${B2_OPTIONS}
	INSTALL_COMMAND ""
)
endif()


ExternalProject_Add(spdlog
	URL ${SPDLOG_SOURCE_ARCHIVE}
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_C_FLAGS=-O2
)


ExternalProject_Add(googletest
	URL ${GOOGLETEST_SOURCE_ARCHIVE}
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_C_FLAGS=-O2
)


# Pass optional -DGIT_COMMIT_DATE=... / -DGIT_COMMIT_SHA=... version overrides through to the ehunter build
set(gitVersionOverrideArgs "")
if(DEFINED GIT_COMMIT_DATE)
	list(APPEND gitVersionOverrideArgs -DGIT_COMMIT_DATE=${GIT_COMMIT_DATE})
endif()
if(DEFINED GIT_COMMIT_SHA)
	list(APPEND gitVersionOverrideArgs -DGIT_COMMIT_SHA=${GIT_COMMIT_SHA})
endif()

ExternalProject_Add(ehunter
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/ehunter
	BUILD_ALWAYS YES
	TEST_COMMAND "ctest"
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_PREFIX_PATH:PATH=${installDir}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		-DCMAKE_CXX_FLAGS=-O2
		${gitVersionOverrideArgs}
		# Passed even when OFF, so that switching an option off also reaches the ehunter build's cache
		-DUSE_SYSTEM_HTSLIB=${USE_SYSTEM_HTSLIB}
		-DSYSTEM_HTSLIB_PREFIX:PATH=${SYSTEM_HTSLIB_PREFIX}
		-DLINK_SYSTEM_HTSLIB_STATICALLY=${LINK_SYSTEM_HTSLIB_STATICALLY}
		-DUSE_SYSTEM_BOOST=${USE_SYSTEM_BOOST}
		-DSYSTEM_BOOST_PREFIX:PATH=${SYSTEM_BOOST_PREFIX}
		${CMAKE_SOURCE_DIR}/ehunter
)


ExternalProject_Add_StepDependencies(ehunter configure
	${ehunterDependencies})

