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


ExternalProject_Add(libdeflate
	URL https://github.com/ebiggers/libdeflate/releases/download/v1.23/libdeflate-1.23.tar.gz
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
	-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	-DCMAKE_INSTALL_PREFIX=${installDir}
	-DLIBDEFLATE_BUILD_SHARED_LIB=OFF
	-DLIBDEFLATE_BUILD_STATIC_LIB=ON
	-DLIBDEFLATE_BUILD_GZIP=OFF
	-DCMAKE_POSITION_INDEPENDENT_CODE=ON
)

ExternalProject_Add(htslib
	BUILD_IN_SOURCE YES
	DEPENDS libdeflate
	URL https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CONFIGURE_COMMAND env CFLAGS=-O2\ -g\ -I${installDir}/include LDFLAGS=-L${installDir}/${CMAKE_INSTALL_LIBDIR} ./configure --prefix=${installDir} --enable-libcurl --enable-s3 --enable-gcs --with-libdeflate
	BUILD_COMMAND $(MAKE) lib-static
)


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
	URL https://downloads.sourceforge.net/project/boost/boost/1.87.0/boost_1_87_0.tar.gz
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	PATCH_COMMAND ${BOOST_PATCH_COMMAND}
	CONFIGURE_COMMAND ./bootstrap.sh --prefix=${installDir}
	BUILD_COMMAND ./b2 install -j8 cxxstd=17 -sNO_ZSTD=1 --with-filesystem --with-system --with-iostreams --with-program_options link=static ${B2_OPTIONS}
	INSTALL_COMMAND ""
)


ExternalProject_Add(spdlog
	URL https://github.com/gabime/spdlog/archive/refs/tags/v1.11.0.tar.gz
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)


ExternalProject_Add(googletest
	URL https://github.com/google/googletest/archive/refs/tags/v1.16.0.tar.gz
	DOWNLOAD_EXTRACT_TIMESTAMP OLD
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
)


ExternalProject_Add(ehunter
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/ehunter
	BUILD_ALWAYS YES
	TEST_COMMAND "ctest"
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${installDir}
		-DCMAKE_PREFIX_PATH:PATH=${installDir}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${CMAKE_SOURCE_DIR}/ehunter
)


ExternalProject_Add_StepDependencies(ehunter configure
	Boost spdlog htslib googletest)

