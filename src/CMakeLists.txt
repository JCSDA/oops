# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_src_files
  oops/oopsNamespaceDoc.h
)

ecbuild_add_library( TARGET     ${PROJECT_NAME}
                     SOURCES    ${oops_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
                     LINKER_LANGUAGE CXX )

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/oops/util/linkedList_i.f
               ${CMAKE_CURRENT_SOURCE_DIR}/oops/util/linkedList_c.f
               DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oops/util )

target_compile_features( ${PROJECT_NAME} PUBLIC cxx_std_14 )

if(OpenMP_FOUND)
    target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX )
endif()
target_link_libraries( ${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran ) # For utils/netcdf_utils_mod
target_link_libraries( ${PROJECT_NAME} PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} )
target_link_libraries( ${PROJECT_NAME} PUBLIC Eigen3::Eigen )
target_link_libraries( ${PROJECT_NAME} PUBLIC eckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC fckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC atlas_f )
target_link_libraries( ${PROJECT_NAME} PUBLIC Boost::boost )
if ( GPTL_FOUND )
	target_link_libraries( ${PROJECT_NAME} PUBLIC GPTL::GPTL )
endif()

# A convenience function to make adding source files from subdirectories easier
# (to keep from having to use the full path every time)
function(oops_target_sources)
  foreach( src IN LISTS ARGN)
    target_sources( ${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${src})
  endforeach()
endfunction()

add_subdirectory( oops )
add_subdirectory( test )

# Stack traces on floating point exceptions
include( backtrace_deps )

target_link_libraries( ${PROJECT_NAME} PUBLIC ${OOPS_STACKTRACE_${OOPS_STACKTRACE_PROVIDER}_LIBS} )
target_compile_definitions( ${PROJECT_NAME} PRIVATE ${OOPS_STACKTRACE_${OOPS_STACKTRACE_PROVIDER}_DEFS} )

# Configure include directory layout for build-tree to match install-tree
set(OOPS_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
target_include_directories(${PROJECT_NAME} INTERFACE
                                $<BUILD_INTERFACE:${OOPS_BUILD_DIR_INCLUDE_PATH}>
                                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                $<INSTALL_INTERFACE:include/oops>)

# Fortran module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                                $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                                $<INSTALL_INTERFACE:${MODULE_DIR}>)

if (${nlohmann_json_FOUND} AND ${nlohmann_json_schema_validator_FOUND})
  message(STATUS "[OOPS] JSON Schema Validator activated.")
  target_link_libraries( oops PRIVATE nlohmann_json nlohmann_json_schema_validator )
  target_compile_definitions( oops PRIVATE OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR )
endif()
