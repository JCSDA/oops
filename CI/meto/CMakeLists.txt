cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(oops-ci VERSION 1.0.0 LANGUAGES C CXX Fortran)

set(ENABLE_MPI ON CACHE BOOL "Compile with MPI")
set(ENABLE_OMP ON CACHE BOOL "Compile with OpenMP")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/jedicmake")
if(NOT DEFINED jedicmake_DIR)
  set(jedicmake_DIR "${CMAKE_BINARY_DIR}/jedicmake")
endif()
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/oops")

include(CTest)
enable_testing()
