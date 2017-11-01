list( APPEND qg_test_input
  testinput/3dvar.json
  testinput/3dfgat.json
  testinput/4densvar.json
  testinput/4dvar.alpha.json
  testinput/4dforcing.json
  testinput/4dhybrid.json
  testinput/4dsaddlepoint.json
  testinput/4dvar.dripcg.json
  testinput/4dvar.drpcglmp.json
  testinput/4dvar.drpfom.json
  testinput/4dvar.drplanczos.json
  testinput/4dvar.ipcg.json
  testinput/4dvar.obsbias.json
  testinput/4dvar.rpcg.json
  testinput/dfi.json
  testinput/forecast.json
  testinput/genenspert.json
  testinput/interfaces.json
  testinput/makeobs3d.json
  testinput/makeobs4d.json
  testinput/makeobs4d12h.json
  testinput/makeobsbias.json
  testinput/test.json
  testinput/truth.json
  testinput/test_op_obs.json
  compare.sh
)

list( APPEND qg_testoutput
  testoutput/3dvar.test
  testoutput/3dfgat.test
  testoutput/4densvar.test
  testoutput/4dvar.alpha.test
  testoutput/4dforcing.test
  testoutput/4dhybrid.test
  testoutput/4dsaddlepoint.test
  testoutput/4dvar.dripcg.test
  testoutput/4dvar.drpcglmp.test
  testoutput/4dvar.drpfom.test
  testoutput/4dvar.drplanczos.test
  testoutput/4dvar.ipcg.test
  testoutput/4dvar.obsbias.test
  testoutput/4dvar.rpcg.test
  testoutput/dfi.test
  testoutput/forecast.test
  testoutput/genenspert.test
  testoutput/makeobs3d.test
  testoutput/makeobs4d.test
  testoutput/makeobs4d12h.test
  testoutput/makeobsbias.test
  testoutput/truth.test
)

# This line copies all files to binary dir
#file( COPY ${qg_test_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${qg_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${qg_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   qg_test_scripts
                       SOURCES_PACK
                       ${qg_test_input}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

ecbuild_add_test( TARGET test_qg_truth
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/truth.json"
                       testoutput/truth.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_forecast
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/forecast.json"
                       testoutput/forecast.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_makeobs3d
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobs3d.json"
                       testoutput/makeobs3d.test
                  DEPENDS qg_makeobs.x )

ecbuild_add_test( TARGET test_qg_makeobs4d
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobs4d.json"
                       testoutput/makeobs4d.test
                  DEPENDS qg_makeobs.x )

ecbuild_add_test( TARGET test_qg_makeobs4d12h
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobs4d12h.json"
                       testoutput/makeobs4d12h.test
                  DEPENDS qg_makeobs.x )

ecbuild_add_test( TARGET test_qg_makeobsbias
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobsbias.json"
                       testoutput/makeobsbias.test
                  DEPENDS qg_makeobs.x )

# Test QG internal classes

# None yet

# Test interface classes with QG

ecbuild_add_test( TARGET  test_qg_geometry
                  BOOST
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_locations
                  BOOST
                  SOURCES executables/TestLocations.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_variables
                  BOOST
                  SOURCES executables/TestVariables.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_geovals
                  BOOST
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_modelaux
                  BOOST
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

#ecbuild_add_test( TARGET  test_qg_modelauxincrement
#                  BOOST
#                  SOURCES executables/TestModelAuxIncrement.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_modelauxcovariance
                  BOOST
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_linearmodel
                  BOOST
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_observationspace
                  BOOST
                  SOURCES executables/TestObservationSpace.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_obsvector
                  BOOST
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsoperator
                  BOOST
                  SOURCES executables/TestObsOperator.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_linobsoperator
                  BOOST
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_obserrorcovar
                  BOOST
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsaux
                  BOOST
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsauxincrement
                  BOOST
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsauxcovariance
                  BOOST
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_localization
                  BOOST
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    qg )

# Test oops classes with QG

#ecbuild_add_test( TARGET test_qg_op_obs.x
#                  BOOST
#                  SOURCES executables/TestSuiteOpObs.cc
#                  ARGS "testinput/test_op_obs.json"
#                  LIBS    qg )

# Test oops applications with QG

ecbuild_add_test( TARGET test_qg_3dvar
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar.json"
                       testoutput/3dvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dfgat
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dfgat.json"
                       testoutput/3dfgat.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drpcglmp
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.drpcglmp.json"
                       testoutput/4dvar.drpcglmp.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drpfom
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.drpfom.json"
                       testoutput/4dvar.drpfom.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.drplanczos.json"
                       testoutput/4dvar.drplanczos.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_dripcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.dripcg.json"
                       testoutput/4dvar.dripcg.test
                  DEPENDS qg_4dvar.x )

if( ${OOPS_SENSITIVE_TESTS} )
ecbuild_add_test( TARGET test_qg_4dvar_ipcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.ipcg.json"
                       testoutput/4dvar.ipcg.test
                  DEPENDS qg_4dvar.x )
endif()

ecbuild_add_test( TARGET test_qg_4dvar_rpcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.rpcg.json"
                       testoutput/4dvar.rpcg.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_obsbias
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.obsbias.json"
                       testoutput/4dvar.obsbias.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_forcing
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dforcing.json"
                       testoutput/4dforcing.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_saddlepoint
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dsaddlepoint.json"
                       testoutput/4dsaddlepoint.test
                  DEPENDS qg_4dvar.x )

if( ${OOPS_SENSITIVE_TESTS} )
ecbuild_add_test( TARGET test_qg_genenspert
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_genpert.x testinput/genenspert.json"
                       testoutput/genenspert.test
                  DEPENDS qg_genpert.x )

ecbuild_add_test( TARGET test_qg_4dvar_hybrid
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dhybrid.json"
                       testoutput/4dhybrid.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4densvar
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4densvar.json"
                       testoutput/4densvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_alpha
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.alpha.json"
                       testoutput/4dvar.alpha.test
                  DEPENDS qg_4dvar.x )
endif()

ecbuild_add_test( TARGET test_qg_dfi
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dfi.x testinput/dfi.json"
                       testoutput/dfi.test
                  DEPENDS qg_dfi.x )

