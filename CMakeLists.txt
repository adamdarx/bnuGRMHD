cmake_minimum_required(VERSION 3.16)

# Project name and source file language
project(GRMHD
        LANGUAGES C CXX)

add_executable(GRMHD main.cpp Metric.cpp)
include_directories(~/eigen)			# 导入Eigen库，默认在主文件夹下


# To use a pre-installed AMReX build, run:
#    cmake -DAMReX_ROOT=/path/to/installdir
# Otherwise cmake will download AMReX from GitHub

if(NOT DEFINED AMReX_ROOT)
  message("-- Download and configure AMReX from GitHub")

  #Download AMReX from GitHub
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)

  FetchContent_Declare(
    amrex_code
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git/
    GIT_TAG        origin/development
    )

  FetchContent_Populate(amrex_code)

  # CMake will read the files in these directories and configure, build
  # and install AMReX.
  add_subdirectory(${amrex_code_SOURCE_DIR} ${amrex_code_BINARY_DIR})

else()

  # Add AMReX
  message("-- Searching for AMReX install directory at ${AMReX_ROOT}")
  find_package(AMReX PATHS ${AMReX_ROOT}/lib/cmake/AMReX/AMReXConfig.cmake)

endif()

target_link_libraries(GRMHD PRIVATE AMReX::amrex)
