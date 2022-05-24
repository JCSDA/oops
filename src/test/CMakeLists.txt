# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND oops_test_src_files
  TestEnvironment.h
  TestFixture.h
)

# oops test input files
list( APPEND oops_test_input
  testinput/random.yaml
  testinput/pushstringvector.yaml
  testinput/atlas_interpolation_interface.yaml
  testinput/unstructured_interpolation.yaml
  testinput/variables.yaml
  testinput/parameters.yaml
  testinput/parameters_older_eckit.yaml
  testinput/empty.yaml
  testinput/mpi.yaml
  testinput/spectrallmp.yaml
  testinput/fft_multiple.yaml
  testinput/hello.yaml
)
# oops test output files
list( APPEND oops_test_output
  testoutput/fft_multiple.test
)

# Macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${FILENAME} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${oops_test_input} )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${oops_test_output} )

oops_target_sources( ${oops_test_src_files} )
install( DIRECTORY test DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oops FILES_MATCHING PATTERN "*.h")

# Configure include directory layout for build-tree to match install-tree
set(OOPS_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(oops_test_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory ${OOPS_BUILD_DIR_INCLUDE_PATH}/oops
                                           COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${OOPS_BUILD_DIR_INCLUDE_PATH}/oops/test)
add_dependencies(oops oops_test_headers)


add_subdirectory( base )
add_subdirectory( assimilation )
add_subdirectory( generic )
add_subdirectory( interface )
add_subdirectory( mpi )
add_subdirectory( util )


ecbuild_add_test( TARGET oops_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

ecbuild_add_executable(
  TARGET  test_oops_base_dummy_run
  SOURCES base/DummyRun.cc
  LIBS oops
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_one
  COMMAND test_oops_base_dummy_run
  ARGS    "testinput/hello.yaml"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_one PROPERTIES
  PASS_REGULAR_EXPRESSION "hello earth"
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_validate_zero
  COMMAND test_oops_base_dummy_run
  ARGS    "--validate-only" "testinput/hello.yaml"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_validate_zero PROPERTIES
  PASS_REGULAR_EXPRESSION "Configuration OK"
)

if (${nlohmann_json_FOUND} AND ${nlohmann_json_schema_validator_FOUND})
  ecbuild_add_test(
    TARGET  test_oops_base_dummy_run_validate_one
    COMMAND test_oops_base_dummy_run
    ARGS    "--validate-only" "testinput/empty.yaml"
    DEPENDS test_oops_base_dummy_run
  )
  set_tests_properties(test_oops_base_dummy_run_validate_one PROPERTIES
    WILL_FAIL true
  )
endif()

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_zero
  COMMAND test_oops_base_dummy_run
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_zero PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_one
  COMMAND test_oops_base_dummy_run
                  ARGS    "--validate-only"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_one PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_two
  COMMAND test_oops_base_dummy_run
  ARGS    "--output-json-schema=x" "testinput/hello.yaml"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_two PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_three
  COMMAND test_oops_base_dummy_run
  ARGS    "--validate-only" "testinput/hello.yaml" "out"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_three PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_help
  COMMAND test_oops_base_dummy_run
  ARGS    "--help"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_help PROPERTIES
  PASS_REGULAR_EXPRESSION "Usages:"
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_h
  COMMAND test_oops_base_dummy_run
  ARGS    "-h"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_h PROPERTIES
  PASS_REGULAR_EXPRESSION "Usages:"
)

ecbuild_add_test( TARGET  test_oops_base_variables
                  SOURCES base/Variables.cc
                  ARGS    "testinput/variables.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_base_posttimer
                  SOURCES base/PostTimer.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops eckit)

ecbuild_add_test( TARGET  test_util_signal_trap
                  SOURCES util/signal_trap.cc
                  LIBS    oops)
set_property( TEST test_util_signal_trap PROPERTY WILL_FAIL TRUE )

ecbuild_add_test( TARGET  test_util_random
                  SOURCES base/Random.cc
                  ARGS    "testinput/random.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_pushstringvector
                  SOURCES base/FCString.cc
                  ARGS    "testinput/pushstringvector.yaml"
                  LIBS    oops )

if (eckit_VERSION VERSION_GREATER_EQUAL 1.16)
  ecbuild_add_test( TARGET  test_util_parameters
                    SOURCES base/Parameters.cc
                    ARGS    "testinput/parameters.yaml"
                    LIBS    oops )
else()
  # Due to changes in eckit 1.16 a different input file is required.
  ecbuild_add_test( TARGET  test_util_parameters
                    SOURCES base/Parameters.cc
                    ARGS    "testinput/parameters_older_eckit.yaml"
                    LIBS    oops )
endif()

ecbuild_add_test( TARGET  test_generic_gc99
                  SOURCES generic/gc99.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_generic_soar
                  SOURCES generic/soar.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_isanypointinvolumeinterior
                  SOURCES util/IsAnyPointInVolumeInterior.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_partialdatetime
                  SOURCES util/PartialDateTime.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_datetime
                  SOURCES util/DateTime.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_duration
                  SOURCES util/Duration.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_intset_parser
                  SOURCES util/IntSetParser.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_scalarormap
                  SOURCES util/ScalarOrMap.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_floatcompare
                  SOURCES util/FloatCompare.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_compositepath
                  SOURCES util/CompositePath.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_stringfunctions
                  SOURCES util/stringFunctions.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_testreference
                  SOURCES util/TestReference.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_mpi_mpi
                  MPI     4
                  SOURCES mpi/mpi.cc
                  ARGS    "testinput/mpi.yaml"
                  LIBS    oops )

#ecbuild_add_test( TARGET  test_oops_unstructured_interpolation
#                  SOURCES generic/InterpolationInterface.cc
#                  ARGS    "testinput/unstructured_interpolation.yaml"
#                  LIBS    oops )

ecbuild_add_test( TARGET  test_fft_multiple
                  SOURCES generic/fft_multiple.cc
                  ARGS    "testinput/fft_multiple.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_algorithms
                  SOURCES util/algorithms.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_comparenvectors
                  SOURCES util/CompareNVectors.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_missingvalues
                  SOURCES util/MissingValues.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_associativecontainers
                  SOURCES util/AssociativeContainers.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_propertiesofnvectors
                  SOURCES util/PropertiesOfNVectors.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_localenvironment
                  SOURCES util/LocalEnvironment.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_typetraits
                  SOURCES util/TypeTraits.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_wildcard
                  SOURCES util/wildcard.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_fullgmres
                  SOURCES assimilation/FullGMRES.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_rotmat
                  SOURCES assimilation/rotmat.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_solvematrixequation
                  SOURCES assimilation/SolveMatrixEquation.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_spectrallmp
                  SOURCES assimilation/SpectralLMP.cc
                  ARGS    "testinput/spectrallmp.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_testvector3d
                  SOURCES assimilation/TestVector3D.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_tridiagsolve
                  SOURCES assimilation/TriDiagSolve.cc
                  ARGS    "testinput/empty.yaml"
                  LIBS    oops )

if( HAVE_FCTEST )

add_fctest( TARGET  test_util_datetime_intfb.x
            SOURCES util/datetime.F90
            LINKER_LANGUAGE Fortran
            CONDITION HAVE_FCTEST
            LIBS    oops )

add_fctest( TARGET  test_util_duration_intfb.x
            SOURCES util/duration.F90
            LINKER_LANGUAGE Fortran
            CONDITION HAVE_FCTEST
            LIBS    oops )

add_fctest( TARGET  test_util_missing_values_f.x
            SOURCES util/missing_values.F90
            LINKER_LANGUAGE Fortran
            CONDITION HAVE_FCTEST
            LIBS    oops )

endif()
