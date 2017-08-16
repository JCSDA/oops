################################################################################
# cmake options:
#
#       -DCMAKE_BUILD_TYPE=Debug|RelWithDebInfo|Release|Production
#       -DCMAKE_INSTALL_PREFIX=/path/to/install
#
#       -DCMAKE_C_COMPILER=gcc
#       -DCMAKE_CXX_COMPILER=g++
#
#       -DCMAKE_PREFIX_PATH=/path/to/jasper:/path/to/any/package/out/of/place
#       -DBUILD_SHARED_LIBS=OFF
################################################################################

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( oops C CXX Fortran)

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake")

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Do compile with MPI" )
set( FIND_OMP ON CACHE BOOL "Do compile with OpenMP" )
set( HAVE_OMP ON ) # Why?
set( NETCDF_F90 ON CACHE BOOL "Do compile with Fortran NetCDF" )

include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_requires_macro_version( 1.6 )

################################################################################
# OOPS project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

# eckit
ecbuild_use_package( PROJECT eckit VERSION 0.12.0 REQUIRED )
include_directories( ${ECKIT_INCLUDE_DIRS} )

# fckit (optional for now)
ecbuild_use_package( PROJECT fckit )
ecbuild_add_option( FEATURE FCTEST
                    DEFAULT ${ENABLE_TESTS}
                    DESCRIPTION "Fortran Unit Testing Framework"
                    REQUIRED_PACKAGES "PROJECT fckit" )
include_directories( ${FCKIT_INCLUDE_DIRS} )

list( APPEND OOPS_TPLS eckit )

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
if( NOT ( BUILD_SITE MATCHES "cc[abt].*" ) )  # Should depend on toolchain, not machine
  ecbuild_find_package( NAME LAPACK REQUIRED )
  include_directories( ${LAPACK_INCLUDE_DIR} )
endif()

# Xml2
find_package( LibXml2 REQUIRED )
include_directories( ${LIBXML2_INCLUDE_DIR} )

# NetCDF
find_package( NetCDF REQUIRED )
include_directories( ${NETCDF_INCLUDE_DIR} )

# OpenMP
find_package( OpenMP REQUIRED )

################################################################################
# Export package info
################################################################################

list( APPEND OOPS_TPLS LAPACK LibXml2 NetCDF OpenMP )
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

ecbuild_add_cxx11_flags()
ecbuild_add_option( FEATURE MPI DEFAULT OFF )

if ( HAVE_MPI AND NOT ECKIT_HAVE_MPI )
  ecbuild_error("OOPS was requested to build with MPI but eckit wasn't built with MPI: ECKIT_HAVE_MPI = ${ECKIT_HAVE_MPI}")
endif()

include( oops_compiler_flags )
include_directories( ${OOPS_INCLUDE_DIRS} ${OOPS_EXTRA_INCLUDE_DIRS} )
add_subdirectory( src )

# Build plugins: The idea is that these become separate builds,
# which find OOPS as an installed package or subproject

add_subdirectory( l95 )
add_subdirectory( qg )

add_subdirectory( ModelX )

ecbuild_dont_pack ( FILES .pydevproject .project .cproject oops.sublime-project )
ecbuild_dont_pack ( FILES .settings tools l95/tools qg/scripts )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME oops )

# print the summary of the configuration
ecbuild_print_summary()
