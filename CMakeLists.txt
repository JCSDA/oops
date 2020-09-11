# (C) Copyright 2009-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

################################################################################
# OOPS
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( oops LANGUAGES CXX Fortran )

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
option( ENABLE_LORENZ95_MODEL "Build LORENZ95 toy model" ON )
option( ENABLE_QG_MODEL "Build QG toy model" ON )
option( ENABLE_OOPS_DOC "Build OOPS documentation" OFF )
option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
option( ENABLE_GPTL "Use GPTL profiling library (if available)" OFF )
option( ENABLE_AUTOPROFILING "Enable function-based autoprofiling with GPTL (if available)" OFF )

include( ${PROJECT_NAME}_compiler_flags )
include( GNUInstallDirs )

## Dependencies
if( ENABLE_MKL )
    find_package( MKL )
endif()
if( MKL_FOUND )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else()
    find_package( LAPACK REQUIRED )
endif()

find_package( Eigen3 REQUIRED NO_MODULE HINTS
    $ENV{Eigen3_ROOT} $ENV{EIGEN3_ROOT} $ENV{Eigen_ROOT} $ENV{EIGEN_ROOT}
    $ENV{Eigen3_PATH} $ENV{EIGEN3_PATH} $ENV{Eigen_PATH} $ENV{EIGEN_PATH})

find_package( OpenMP COMPONENTS CXX Fortran )
find_package( MPI REQUIRED COMPONENTS CXX Fortran )
find_package( NetCDF REQUIRED COMPONENTS Fortran )
find_package( Boost 1.64.0 REQUIRED)
find_package( eckit 1.4.0 REQUIRED COMPONENTS MPI )
find_package( fckit 0.6.3 REQUIRED )
if(OpenMP_FOUND)
    find_package( atlas 0.18.1 REQUIRED COMPONENTS OMP OMP_Fortran )
else()
    find_package( atlas 0.18.1 REQUIRED )
endif()

#Optional Dependencies
find_package( jedicmake QUIET )
if( jedicmake_FOUND AND ENABLE_GPTL )
    find_package( GPTL )
    add_compile_definitions(ENABLE_GPTL)
endif()

if ( GPTL_FOUND AND ENABLE_AUTOPROFILING )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_AUTOPROFILING_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_AUTOPROFILING_FLAGS}")
endif()

## Sources
add_subdirectory( src )
add_subdirectory( tools )

if( ENABLE_OOPS_DOC )
    add_subdirectory( Documents )
endif()

if( ENABLE_LORENZ95_MODEL )
  add_subdirectory( l95 )
endif()
if( ENABLE_QG_MODEL )
  add_subdirectory( qg )
endif()

## Package Config

ecbuild_install_project( NAME ${PROJECT_NAME} )

#Export Targets
set(EXPORTED_TARGET_NAMESPACE "") #TODO: set to "oops::" once all packages are ready for namespaced targets
set(INSTALL_TARGETS_PATH lib/cmake/${PROJECT_NAME})
get_filename_component(TARGETS_FILE "${PROJECT_TARGETS_FILE}" NAME)
# export(EXPORT ${PROJECT_NAME}-targets NAMESPACE "${EXPORTED_TARGET_NAMESPACE}" FILE ${PROJECT_TARGETS_FILE}) #Build-tree
#install(EXPORT ${PROJECT_NAME}-targets NAMESPACE "${EXPORTED_TARGET_NAMESPACE}" DESTINATION ${INSTALL_TARGETS_PATH}) #Install-tree

#Export package config for build and install trees
set( BINDIR_BUILD_EXPORT ${CMAKE_BINARY_DIR}/bin ) #Exported binary (tools) directory locations for build-tree
file( MAKE_DIRECTORY ${BINDIR_BUILD_EXPORT} )
set( BINDIR_INSTALL_EXPORT ${CMAKE_INSTALL_FULL_BINDIR} ) #Exported binary (tools) directory locations for install-tree
ecbuild_generate_project_config(${PROJECT_NAME}-config.cmake.in PATH_VARS BINDIR_BUILD_EXPORT BINDIR_INSTALL_EXPORT )

ecbuild_print_summary()
