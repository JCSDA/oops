
# (C) Copyright 2017 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# OOPS bundle
#

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( oops-bundle VERSION 1.0.0 LANGUAGES C CXX Fortran )

find_package(ecbuild REQUIRED)
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( ecbuild_bundle )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()

ecbuild_bundle( PROJECT oops  GIT "https://github.com/jcsda-internal/oops.git" )

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
