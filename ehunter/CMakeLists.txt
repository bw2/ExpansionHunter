cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter LANGUAGES CXX)

# Git commit date (for the version string) and short git commit SHA (reported in the output JSON as the
# "Version" field). Each is taken from the first of these that provides it:
#   1. -DGIT_COMMIT_DATE=... / -DGIT_COMMIT_SHA=... passed to cmake.
#   2. GitArchiveVersion.txt, which `git archive` fills in (export-subst in .gitattributes), so a GitHub
#      archive/<sha>.tar.gz download, which has no .git, still gets the real values. In a git checkout the
#      file holds the unexpanded $Format:...$ placeholders, which don't match the regexes below.
#   3. git, when building from a git checkout.
#   4. "unknown".
file(STRINGS ${CMAKE_SOURCE_DIR}/GitArchiveVersion.txt GIT_ARCHIVE_VERSION_LINES)

if(NOT DEFINED GIT_COMMIT_DATE)
    if("${GIT_ARCHIVE_VERSION_LINES}" MATCHES "commit_date=([0-9]+-[0-9]+-[0-9]+)")
        set(GIT_COMMIT_DATE "${CMAKE_MATCH_1}")
    else()
        execute_process(
            COMMAND git log -1 --format=%cs
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_DATE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif()
endif()
if(NOT GIT_COMMIT_DATE)
    set(GIT_COMMIT_DATE "unknown")
endif()
add_compile_definitions(GIT_COMMIT_DATE="${GIT_COMMIT_DATE}")

if(NOT DEFINED GIT_COMMIT_SHA)
    if("${GIT_ARCHIVE_VERSION_LINES}" MATCHES "commit_sha=([0-9a-f]+)")
        set(GIT_COMMIT_SHA "${CMAKE_MATCH_1}")
    else()
        execute_process(
            COMMAND git rev-parse --short HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_SHA
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif()
endif()
if(NOT GIT_COMMIT_SHA)
    set(GIT_COMMIT_SHA "unknown")
endif()
add_compile_definitions(GIT_COMMIT_SHA="${GIT_COMMIT_SHA}")


set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_BUILD_TYPE Debug)

if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
            "Choose the type of build (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif ()

message(STATUS "BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# The htslib and Boost search results below are cached, so changing these options in an existing build directory
# would silently keep the libraries found before. Require a fresh build directory instead.
set(systemLibrarySettings "USE_SYSTEM_HTSLIB=${USE_SYSTEM_HTSLIB} SYSTEM_HTSLIB_PREFIX=${SYSTEM_HTSLIB_PREFIX} \
LINK_SYSTEM_HTSLIB_STATICALLY=${LINK_SYSTEM_HTSLIB_STATICALLY} USE_SYSTEM_BOOST=${USE_SYSTEM_BOOST} \
SYSTEM_BOOST_PREFIX=${SYSTEM_BOOST_PREFIX}")
if (DEFINED CONFIGURED_SYSTEM_LIBRARY_SETTINGS AND NOT CONFIGURED_SYSTEM_LIBRARY_SETTINGS STREQUAL systemLibrarySettings)
    message(FATAL_ERROR "The installed-library options changed since this build directory was configured "
            "(was: ${CONFIGURED_SYSTEM_LIBRARY_SETTINGS}; now: ${systemLibrarySettings}). Use a fresh build directory.")
endif ()

set(Boost_USE_STATIC_LIBS ON)
if (USE_SYSTEM_BOOST AND SYSTEM_BOOST_PREFIX)
    set(Boost_ROOT ${SYSTEM_BOOST_PREFIX})
endif ()
find_package(Boost 1.84 REQUIRED COMPONENTS program_options filesystem system iostreams)


add_compile_options(-Wall -Werror -Wextra -O2)

enable_testing()

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)
find_package(spdlog REQUIRED)

if (NOT APPLE)
    find_package(OpenSSL REQUIRED)
endif ()

if (USE_SYSTEM_HTSLIB)
    # By default link the installed shared htslib, which brings its own dependencies (libdeflate, libcurl, ...), so
    # libdeflate is not linked here. With LINK_SYSTEM_HTSLIB_STATICALLY, link libhts.a instead. The headers are
    # taken from the same prefix as the library so that the two always match. The installed ExpansionHunter binary
    # gets no stored library search path (CMake's default install behavior, kept on purpose), so a shared htslib
    # in a custom prefix must be on LD_LIBRARY_PATH when running it.
    if (LINK_SYSTEM_HTSLIB_STATICALLY)
        set(htslibLibraryName libhts.a)
    else ()
        set(htslibLibraryName hts)
    endif ()
    if (SYSTEM_HTSLIB_PREFIX)
        find_path(HTSLIB_INCLUDE_DIR htslib/hts.h PATHS ${SYSTEM_HTSLIB_PREFIX}/include NO_DEFAULT_PATH)
        # lib/<architecture> covers Debian/Ubuntu's multiarch layout (e.g. /usr/lib/x86_64-linux-gnu)
        find_library(htslib ${htslibLibraryName}
                PATHS ${SYSTEM_HTSLIB_PREFIX}/lib/${CMAKE_LIBRARY_ARCHITECTURE} ${SYSTEM_HTSLIB_PREFIX}/lib ${SYSTEM_HTSLIB_PREFIX}/lib64
                NO_DEFAULT_PATH)
    else ()
        find_path(HTSLIB_INCLUDE_DIR htslib/hts.h)
        find_library(htslib ${htslibLibraryName})
    endif ()
    if (NOT HTSLIB_INCLUDE_DIR OR NOT htslib)
        message(FATAL_ERROR "USE_SYSTEM_HTSLIB is ON but htslib was not found "
                "(headers: ${HTSLIB_INCLUDE_DIR}, library: ${htslib}). Set SYSTEM_HTSLIB_PREFIX to its install prefix.")
    endif ()

    # HTS_VERSION is defined from htslib 1.10 onward, which is also the first release with sam_index_load3()
    file(STRINGS ${HTSLIB_INCLUDE_DIR}/htslib/hts.h HTS_VERSION_LINE REGEX "^#define HTS_VERSION [0-9]+")
    string(REGEX MATCH "[0-9]+" HTS_VERSION_NUMBER "${HTS_VERSION_LINE}")
    if (NOT HTS_VERSION_NUMBER OR HTS_VERSION_NUMBER LESS 101000)
        message(FATAL_ERROR "htslib 1.10 or newer is required, but ${HTSLIB_INCLUDE_DIR}/htslib/hts.h is older")
    endif ()
    message(STATUS "Using installed htslib: ${htslib} (HTS_VERSION ${HTS_VERSION_NUMBER})")
    set(htslibDependencyLibraries "")

    # A static htslib needs its dependencies linked here. The static_libs line of the htslib.pc installed next to
    # libhts.a lists them (e.g. "-lpthread -lz -lm -lbz2 -llzma -ldeflate -lcurl", plus -lhtscodecs for builds
    # with an external htscodecs, such as Debian's). Those already linked below are skipped, the C runtime's own
    # libraries (m, dl, rt) are linked by name, and each remaining one is looked up, preferring a static library and
    # looking next to htslib first. Without htslib.pc, only libdeflate is added, when it is found.
    if (LINK_SYSTEM_HTSLIB_STATICALLY)
        get_filename_component(htslibLibraryDir ${htslib} DIRECTORY)
        set(htslibPkgConfigFile ${htslibLibraryDir}/pkgconfig/htslib.pc)
        set(htslibLinkedBelow pthread z bz2 lzma curl crypto)
        set(cRuntimeLibraries m dl rt)
        if (EXISTS ${htslibPkgConfigFile})
            file(STRINGS ${htslibPkgConfigFile} htslibStaticLibsLine REGEX "^static_libs=")
            string(REGEX MATCHALL "-l[^ ]+" htslibStaticLibFlags "${htslibStaticLibsLine}")
            string(REPLACE "-l" "" htslibStaticLibNames "${htslibStaticLibFlags}")
            list(REMOVE_ITEM htslibStaticLibNames ${htslibLinkedBelow})
            set(htslibStaticLibsRequired YES)
        else ()
            set(htslibStaticLibNames deflate)
            set(htslibStaticLibsRequired NO)
        endif ()
        foreach (libraryName ${htslibStaticLibNames})
            if (libraryName IN_LIST cRuntimeLibraries)
                list(APPEND htslibDependencyLibraries ${libraryName})
                continue()
            endif ()
            find_library(HTSLIB_DEPENDENCY_${libraryName} NAMES lib${libraryName}.a ${libraryName} HINTS ${htslibLibraryDir})
            if (HTSLIB_DEPENDENCY_${libraryName})
                message(STATUS "Linking ${libraryName} for the static htslib: ${HTSLIB_DEPENDENCY_${libraryName}}")
                list(APPEND htslibDependencyLibraries ${HTSLIB_DEPENDENCY_${libraryName}})
            elseif (htslibStaticLibsRequired)
                message(FATAL_ERROR "The static htslib needs lib${libraryName} (listed in ${htslibPkgConfigFile}), but it was not found")
            else ()
                message(STATUS "${libraryName} not found; linking the static htslib only works if it was built without it")
            endif ()
        endforeach ()
    endif ()
else ()
    find_package(libdeflate REQUIRED)
    find_path(HTSLIB_INCLUDE_DIR htslib/hts.h)
    find_library(htslib libhts.a)
    find_library(htslib hts)
    set(htslibDependencyLibraries libdeflate::libdeflate_static)
endif ()
# Recorded only once the htslib and Boost checks above have passed, so that fixing a failed configure (e.g. setting
# SYSTEM_HTSLIB_PREFIX after htslib was not found) does not trip the changed-settings check
set(CONFIGURED_SYSTEM_LIBRARY_SETTINGS "${systemLibrarySettings}" CACHE INTERNAL "")

add_subdirectory(thirdparty/graph-tools-master-0cd9399)
# graph-tools includes spdlog headers without linking spdlog, so pass it spdlog's include directory explicitly
# (otherwise it only finds them when spdlog shares an include directory with Boost)
target_include_directories(graphtools SYSTEM PRIVATE $<TARGET_PROPERTY:spdlog::spdlog,INTERFACE_INCLUDE_DIRECTORIES>)

# Embed the genotype-quality model file (GQ_MODEL_FILE, set below) into the
# binary as a generated byte-array source. Regenerated whenever the model file or the
# generator script changes.
set(GQ_MODEL_FILE ${CMAKE_CURRENT_SOURCE_DIR}/data/genotype_quality_model_from_HG002_and_CHM1_CHM13_plus50.20261008.json.gz)
# Basename of the embedded model file (without .json/.json.gz), reported in the output JSON as
# GenotypeQualityModelVersion so calls can be traced to the model that produced them.
get_filename_component(GQ_MODEL_NAME ${GQ_MODEL_FILE} NAME)
string(REGEX REPLACE "\\.json(\\.gz)?$" "" GQ_MODEL_NAME ${GQ_MODEL_NAME})
add_compile_definitions(GQ_MODEL_NAME="${GQ_MODEL_NAME}")
set(GQ_MODEL_GENERATED ${CMAKE_CURRENT_BINARY_DIR}/GenotypeQualityModelData.cpp)
add_custom_command(
        OUTPUT ${GQ_MODEL_GENERATED}
        COMMAND ${CMAKE_COMMAND} -DIN_FILE=${GQ_MODEL_FILE} -DOUT_FILE=${GQ_MODEL_GENERATED}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedModel.cmake
        DEPENDS ${GQ_MODEL_FILE} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedModel.cmake
        COMMENT "Embedding genotype-quality model from ${GQ_MODEL_FILE}"
        VERBATIM)

add_library(ExpansionHunterLib
        ${GQ_MODEL_GENERATED}
        alignment/AlignmentClassifier.hh alignment/AlignmentClassifier.cpp
        alignment/AlignmentFilters.hh alignment/AlignmentFilters.cpp
        alignment/ClassifierOfAlignmentsToVariant.hh alignment/ClassifierOfAlignmentsToVariant.cpp
        alignment/GraphVariantAlignmentStats.hh alignment/GraphVariantAlignmentStats.cpp
        alignment/GreedyAlignmentIntersector.hh alignment/GreedyAlignmentIntersector.cpp
        alignment/HighQualityBaseRunFinder.hh alignment/HighQualityBaseRunFinder.cpp
        alignment/OperationsOnAlignments.hh alignment/OperationsOnAlignments.cpp
        alignment/OrientationPredictor.hh alignment/OrientationPredictor.cpp
        alignment/SoftclippingAligner.hh alignment/SoftclippingAligner.cpp
        core/Common.hh core/Common.cpp
        core/ConcurrentQueue.hh
        core/CountTable.hh core/CountTable.cpp
        core/GenomicRegion.hh core/GenomicRegion.cpp
        core/HtsHelpers.hh core/HtsHelpers.cpp
        core/Parameters.hh core/Parameters.cpp
        core/Reference.hh core/Reference.cpp
        core/ReferenceContigInfo.hh core/ReferenceContigInfo.cpp
        core/LocusStats.hh core/LocusStats.cpp
        core/LogSum.hh
        core/Read.hh core/Read.cpp
        core/ReadPairs.hh core/ReadPairs.cpp
        core/ReadSupportCalculator.hh core/ReadSupportCalculator.cpp
        core/ThreadPool.hh
        core/WeightedPurityCalculator.hh core/WeightedPurityCalculator.cpp
        core/RepeatPurity.hh core/RepeatPurity.cpp
        genotyping/AlignMatrix.hh genotyping/AlignMatrix.cpp
        genotyping/AlignMatrixFiltering.hh genotyping/AlignMatrixFiltering.cpp
        genotyping/AlleleChecker.hh genotyping/AlleleChecker.cpp
        genotyping/FragLogliks.hh genotyping/FragLogliks.cpp
        genotyping/OneAlleleStrGenotyper.hh genotyping/OneAlleleStrGenotyper.cpp
        genotyping/RepeatGenotype.hh genotyping/RepeatGenotype.cpp
        genotyping/SmallVariantGenotype.hh genotyping/SmallVariantGenotype.cpp
        genotyping/SmallVariantGenotyper.hh genotyping/SmallVariantGenotyper.cpp
        genotyping/StrAlign.hh genotyping/StrAlign.cpp
        genotyping/StrGenotyper.hh genotyping/StrGenotyper.cpp
        genotyping/TwoAlleleStrGenotyper.hh genotyping/TwoAlleleStrGenotyper.cpp
        genotype_quality/GenotypeQualityModel.hh genotype_quality/GenotypeQualityModel.cpp
        genotype_quality/GenotypeQualityFeatures.hh genotype_quality/GenotypeQualityFeatures.cpp
        genotype_quality/GenotypeQualityAnnotator.hh genotype_quality/GenotypeQualityAnnotator.cpp
        io/BamletWriter.hh io/BamletWriter.cpp
        io/CatalogLoading.hh io/CatalogLoading.cpp
        io/GraphBlueprint.hh io/GraphBlueprint.cpp
        io/JsonWriter.hh io/JsonWriter.cpp
        io/IterativeJsonWriter.hh io/IterativeJsonWriter.cpp
        io/IterativeVcfWriter.hh io/IterativeVcfWriter.cpp
        io/StreamingOutputMerge.hh io/StreamingOutputMerge.cpp
        io/ResumeCheckpoint.hh io/ResumeCheckpoint.cpp
        io/LocusSpecDecoding.hh io/LocusSpecDecoding.cpp
        io/ParameterLoading.hh io/ParameterLoading.cpp
        io/RegionGraph.hh io/RegionGraph.cpp
        io/SampleStats.hh io/SampleStats.cpp
        io/StringUtils.hh io/StringUtils.cpp
        io/VcfHeader.hh io/VcfHeader.cpp
        io/VcfWriter.hh io/VcfWriter.cpp
        io/VcfWriterHelpers.hh io/VcfWriterHelpers.cpp
        sample/AnalyzerFinder.hh sample/AnalyzerFinder.cpp
        sample/GenomeMask.hh sample/GenomeMask.cpp
        sample/GenomeQueryCollection.hh sample/GenomeQueryCollection.cpp
        sample/HtsFileSeeker.hh sample/HtsFileSeeker.cpp
        sample/HtsFileStreamer.hh sample/HtsFileStreamer.cpp
        sample/HtsSeekingSampleAnalysis.hh sample/HtsSeekingSampleAnalysis.cpp
        sample/HtsStreamingReadPairQueue.hh sample/HtsStreamingReadPairQueue.cpp
        sample/HtsStreamingSampleAnalysis.hh sample/HtsStreamingSampleAnalysis.cpp
        sample/HtsLowMemStreamingSampleAnalysis.hh sample/HtsLowMemStreamingSampleAnalysis.cpp
        sample/HtsLowMemStreamingHelpers.hh sample/HtsLowMemStreamingHelpers.cpp
        sample/MateExtractor.hh sample/MateExtractor.cpp
        )

set(CTPL_INCLUDE_DIRS thirdparty/ctpl/ctpl-0.0.2)
target_include_directories(ExpansionHunterLib PUBLIC ${CMAKE_SOURCE_DIR})
# spdlog's include directory (the downloaded spdlog and googletest) goes first so that other spdlog or googletest
# headers in a shared prefix such as /opt/homebrew/include can't shadow them
target_include_directories(ExpansionHunterLib SYSTEM PUBLIC
        $<TARGET_PROPERTY:spdlog::spdlog,INTERFACE_INCLUDE_DIRECTORIES>
        ${HTSLIB_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
        ${LIBLZMA_INCLUDE_DIRS}
        ${CURL_INCLUDE_DIRS}
        ${CTPL_INCLUDE_DIRS}
        )

# Set static linking of gcc standard libraries to simplify binary distribution
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(STATIC_FLAGS -static-libgcc -static-libstdc++)
endif ()

target_link_libraries(ExpansionHunterLib PUBLIC
        ${STATIC_FLAGS}
        graphtools
        ${htslib}
        ${Boost_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
        ${CURL_LIBRARIES}
        ${htslibDependencyLibraries}
        ZLIB::ZLIB
        BZip2::BZip2
        Threads::Threads
        spdlog::spdlog)

if (NOT APPLE)
    target_link_libraries(ExpansionHunterLib PUBLIC
            OpenSSL::Crypto)
endif ()

add_executable(ExpansionHunter
        app/ExpansionHunter.cpp
        )
target_link_libraries(ExpansionHunter ExpansionHunterLib)

# Link the mimalloc allocator when available. On allocation-heavy genotyping (long alleles / deep pileups)
# this is ~30% faster and uses less peak RSS than the system allocator, with byte-identical output. Optional
# so builds without mimalloc still succeed (just without the speedup); install mimalloc to enable it.
# Linux only: mimalloc's static malloc override is reliable on Linux but aborts at runtime on macOS, so it is
# disabled on Apple (the production platform is Linux anyway).
if (NOT APPLE)
    find_package(mimalloc QUIET)
    if (TARGET mimalloc-static)
        target_link_libraries(ExpansionHunter mimalloc-static)
        message(STATUS "Linking mimalloc allocator (mimalloc-static)")
    elseif (TARGET mimalloc)
        target_link_libraries(ExpansionHunter mimalloc)
        message(STATUS "Linking mimalloc allocator (mimalloc, shared)")
    else ()
        message(STATUS "mimalloc not found; using the system allocator (install mimalloc for ~30% faster genotyping)")
    endif ()
endif ()

add_executable(UnitTests
        tests/AlignMatrixTest.cpp
        tests/AlignmentBufferTest.cpp
        tests/AnalysisModeTest.cpp
        tests/AlignmentClassifierTest.cpp
        tests/AlignmentSummaryTest.cpp
        tests/AlleleCheckerTest.cpp
        tests/CatalogLoadingKnownMotifsTest.cpp
        tests/CatalogLoadingPlotTest.cpp
        tests/ChromosomeStrideTest.cpp
        tests/ClassifierOfAlignmentsToVariantTest.cpp
        tests/ConcurrentQueueTest.cpp
        tests/ConsensusSequenceTest.cpp
        tests/CountTableTest.cpp
        tests/FragLogliksTest.cpp
        tests/GenomeMaskTest.cpp
        tests/GenotypeQualityModelTest.cpp
        tests/GenotypeQualityFeaturesTest.cpp
        tests/GenotypeQualityAnnotatorTest.cpp
        tests/GenomicRegionTest.cpp
        tests/GraphAlignmentOperationsTest.cpp
        tests/GraphBlueprintTest.cpp
        tests/GreedyAlignmentIntersectorTest.cpp
        tests/HighQualityBaseRunFinderTest.cpp
        tests/HtsLowMemStreamingHelpersTest.cpp
        tests/LocusSpecificationPlotTest.cpp
        tests/LocusStatsTest.cpp
        tests/MateExtractorMergeTest.cpp
        tests/MotifCompositionTest.cpp
        tests/ReadSupportCalculatorTest.cpp
        tests/ReadTest.cpp
        tests/RegionGraphTest.cpp
        tests/RepeatAnalyzerTest.cpp
        tests/RepeatPurityTest.cpp
        tests/ResumeTest.cpp
        tests/RepeatGenotypeTest.cpp
        tests/ReviewerAlignsTest.cpp
        tests/ReviewerFragLenFilterTest.cpp
        tests/ReviewerGenotypePathsTest.cpp
        tests/ReviewerLanePlotTest.cpp
        tests/ReviewerOriginTest.cpp
        tests/ReviewerPhasingTest.cpp
        tests/ReviewerProjectionTest.cpp
        tests/ReviewerWorkflowTest.cpp
        tests/RFC1MotifAnalysisUtilTest.cpp
        tests/SmallVariantGenotyperTest.cpp
        tests/SoftclippingAlignerTest.cpp
        tests/StrAlignTest.cpp
        tests/StrGenotyperTest.cpp
        tests/UnitTests.cpp
        tests/VcfHeaderTest.cpp
        tests/WeightedPurityCalculatorTest.cpp
        )
add_subdirectory(locus)
add_subdirectory(reviewer)

target_link_libraries(UnitTests ExpansionHunterLib GTest::GTest GTest::Main)
target_compile_definitions(UnitTests PRIVATE EH_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")

add_test(NAME UnitTests COMMAND UnitTests)

install(TARGETS ExpansionHunter UnitTests RUNTIME DESTINATION bin)
