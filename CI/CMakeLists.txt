  
# (C) Copyright 2017 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# OOPS bundle
#

project( oops-bundle C CXX Fortran )

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake;${CMAKE_MODULE_PATH}")
include( ecbuild_bundle )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()

ecbuild_requires_macro_version( 2.7 )

ecbuild_add_option( FEATURE OMP
                    DEFAULT ON
                    DESCRIPTION "Use OpenMP" )

ecbuild_bundle( PROJECT fckit GIT "https://github.com/JCSDA/fckit.git" )
ecbuild_bundle( PROJECT atlas GIT "https://github.com/JCSDA/atlas.git" )
ecbuild_bundle( PROJECT oops GIT "https://github.com/JCSDA/oops.git" )

ecbuild_bundle_finalize()
