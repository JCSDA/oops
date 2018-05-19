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

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( oops C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH} )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.7 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()

set( OOPS_LINKER_LANGUAGE CXX )

set( OOPS_SENSITIVE_TESTS true )
if ( CMAKE_CXX_COMPILER_ID MATCHES "Cray" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
  set( OOPS_SENSITIVE_TESTS false )
endif()

################################################################################
# Dependencies
################################################################################

# Boost
set( Boost_MINIMUM_VERSION "1.47" )
find_package( Boost ${Boost_MINIMUM_VERSION} REQUIRED
              COMPONENTS unit_test_framework )
include_directories( ${Boost_INCLUDE_DIR} )

# LAPACK
ecbuild_add_option( FEATURE LAPACK
                    DESCRIPTION "Linear Algebra PACKage"
                    REQUIRED_PACKAGES "LAPACK REQUIRED" )
if (NOT HAVE_LAPACK)
  ecbuild_error("LAPACK_LIBRARIES not defined")
endif()

# NetCDF
set( NETCDF_F90 ON CACHE BOOL "Compile with Fortran NetCDF" )
find_package( NetCDF REQUIRED )
include_directories( ${NETCDF_INCLUDE_DIR} )

# eckit
ecbuild_use_package( PROJECT eckit VERSION 0.18.0 REQUIRED )
include_directories( ${ECKIT_INCLUDE_DIRS} )

# fckit
ecbuild_use_package( PROJECT fckit VERSION 0.4.1 REQUIRED )
include_directories( ${FCKIT_INCLUDE_DIRS} )

################################################################################
# Export package info
################################################################################

list( APPEND OOPS_TPLS LAPACK NetCDF eckit fckit )

set( OOPS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
                       ${CMAKE_Fortran_MODULE_DIRECTORY} )
set( OOPS_LIBRARIES oops )

get_directory_property( OOPS_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${OOPS_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND OOPS_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND OOPS_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND OOPS_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

include( oops_compiler_flags )
include_directories( ${OOPS_INCLUDE_DIRS} ${OOPS_EXTRA_INCLUDE_DIRS} )
add_subdirectory( src )
add_subdirectory( tools )

# Build example applications: these show how to use OOPS but are not part of OOPS

add_subdirectory( l95 )
add_subdirectory( qg )

ecbuild_dont_pack ( FILES .pydevproject .project .cproject oops.sublime-project .settings )

if(ECBUILD_INSTALL_FORTRAN_MODULES)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR} )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME oops )

# print the summary of the configuration
ecbuild_print_summary()
