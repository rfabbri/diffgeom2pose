# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)
set(CMAKE_CXX_STANDARD "11")
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

#------------------------------------------------------------------------------

## TODO: use add_compile_options(-Wall) etc
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall" CACHE STRING "" FORCE)

# Crucial flags for performance
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native" CACHE STRING "" FORCE)

# add_compile_options(
#     "$<$<CONFIG:RELEASE>:-march=native -ffast-math>"
#     "$<$<CONFIG:DEBUG>:-Wall>"
# )

# This one is the most important one - order of magnitude faster
# See this: https://gcc.gnu.org/ml/gcc/2001-07/msg02150.html
# if GCC 8
#   -march=native -freciprocal-math -fno-trapping-math -fno-math-errno -ffinite-math-only -fno-signaling-nans -fno-rounding-math -fcx-limited-range -fexcess-precision=fast
#   problem with gcc 8 is -fno-signed-zeros
# else
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" CACHE STRING "" FORCE)

# --fast-math breaks into these:
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fexcess-precision=fast")

# This is the --fast-math flag that does a big job but sometimes gives wrong
# results in GCC8: -funsafe-math-optimizations
# of the flags implied by this one, the most impact comes from
# -fno-signed-zeros, which is also what breaks GCC8 build
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-signed-zeros -fno-trapping-math -fassociative-math -freciprocal-math")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

#-------------------------------------------------------------------------------

project(DIFFGEOM2POSE #Project name must be all caps to have properly generated PROJECT_VERSION_* variables
    VERSION 1
    DESCRIPTION "Absolute Pose from 2 feature points (substitutes P3P)"
    LANGUAGES CXX C)
include(CMakeDependentOption)
include(${CMAKE_CURRENT_SOURCE_DIR}/config/cmake/DIFFGEOM2POSEInitializeBuildType.cmake)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wno-undefined-var-template HAS_NO_UNDEFINED_VAR_TEMPLATE)
if( HAS_NO_UNDEFINED_VAR_TEMPLATE )
  add_definitions( -Wno-undefined-var-template )
endif()

# Include sub-projects.
add_subdirectory ("cmd")
