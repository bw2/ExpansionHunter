cmake_minimum_required(VERSION 3.13)
project(ExpansionHunter LANGUAGES CXX)

# Get git commit date for version string
execute_process(
    COMMAND git log -1 --format=%cs
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT GIT_COMMIT_DATE)
    set(GIT_COMMIT_DATE "unknown")
endif()
add_compile_definitions(GIT_COMMIT_DATE="${GIT_COMMIT_DATE}")


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

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.84 REQUIRED COMPONENTS program_options filesystem system iostreams)


add_compile_options(-Wall -Werror -Wextra -O2)

enable_testing()

find_package(libdeflate REQUIRED)
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

find_library(htslib libhts.a)
find_library(htslib hts)

add_subdirectory(thirdparty/graph-tools-master-0cd9399)

add_library(ExpansionHunterLib
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
        io/BamletWriter.hh io/BamletWriter.cpp
        io/CatalogLoading.hh io/CatalogLoading.cpp
        io/GraphBlueprint.hh io/GraphBlueprint.cpp
        io/JsonWriter.hh io/JsonWriter.cpp
        io/IterativeJsonWriter.hh io/IterativeJsonWriter.cpp
        io/IterativeVcfWriter.hh io/IterativeVcfWriter.cpp
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
        sample/IndexBasedDepthEstimate.hh sample/IndexBasedDepthEstimate.cpp
        sample/MateExtractor.hh sample/MateExtractor.cpp
        )

set(CTPL_INCLUDE_DIRS thirdparty/ctpl/ctpl-0.0.2)
target_include_directories(ExpansionHunterLib PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(ExpansionHunterLib SYSTEM PUBLIC
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
        libdeflate::libdeflate_static
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

add_executable(UnitTests
        tests/AlignMatrixTest.cpp
        tests/AlignmentBufferTest.cpp
        tests/AlignmentClassifierTest.cpp
        tests/AlignmentSummaryTest.cpp
        tests/AlleleCheckerTest.cpp
        tests/CatalogLoadingPlotTest.cpp
        tests/ClassifierOfAlignmentsToVariantTest.cpp
        tests/ConcurrentQueueTest.cpp
        tests/ConsensusSequenceTest.cpp
        tests/CountTableTest.cpp
        tests/FragLogliksTest.cpp
        tests/GenomeMaskTest.cpp
        tests/GenomicRegionTest.cpp
        tests/GraphAlignmentOperationsTest.cpp
        tests/GraphBlueprintTest.cpp
        tests/GreedyAlignmentIntersectorTest.cpp
        tests/HighQualityBaseRunFinderTest.cpp
        tests/LocusSpecificationPlotTest.cpp
        tests/LocusStatsTest.cpp
        tests/ReadSupportCalculatorTest.cpp
        tests/ReadTest.cpp
        tests/RegionGraphTest.cpp
        tests/RepeatAnalyzerTest.cpp
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
        tests/WeightedPurityCalculatorTest.cpp
        )
add_subdirectory(locus)
add_subdirectory(reviewer)

target_link_libraries(UnitTests ExpansionHunterLib GTest::GTest GTest::Main)

add_test(NAME UnitTests COMMAND UnitTests)

install(TARGETS ExpansionHunter UnitTests RUNTIME DESTINATION bin)
