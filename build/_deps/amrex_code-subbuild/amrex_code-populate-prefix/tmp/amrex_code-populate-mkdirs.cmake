# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src"
  "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build"
  "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix"
  "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/tmp"
  "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp"
  "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src"
  "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-subbuild/amrex_code-populate-prefix/src/amrex_code-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
