list( APPEND l95_test_input
  testinput/3dvar.json
  testinput/3dfgat.json
  testinput/4dvar.alpha.json
  testinput/4dvar.drgmresr.json
  testinput/4dvar.dripcg.json
  testinput/4dvar.dripcgqn.json
  testinput/4dvar.drpcg.json
  testinput/4dvar.drpcgqn.json
  testinput/4dvar.drplanczos.json
  testinput/4dvar.drplanclmp.json
  testinput/4dvar.gmresr.json
  testinput/4dvar.ipcg.json
  testinput/4dvar.pcg.json
  testinput/4dvar.planczos.json
  testinput/4dvar.rpcg.json
  testinput/4dvar.rplanczos.json
  testinput/4dvar.modbias.json
  testinput/4dvar.obsbias.json
  testinput/4dforcing.json
  testinput/4dsaddlepoint.json
  testinput/4dhybrid.json
  testinput/4densvar.json
  testinput/forecast.json
  testinput/genenspert.json
  testinput/hofx.json
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

list( APPEND l95_test_data
  test.an.2010-01-01T00:00:00Z.l95
  truth.an.2010-01-01T00:00:00Z.l95
  truth.fc.2010-01-01T00:00:00Z.PT3H
  truth.2010-01-02T00:00:00Z.obt
)

list( APPEND l95_testoutput
  testoutput/3dvar.test
  testoutput/3dfgat.test
  testoutput/4densvar.test
  testoutput/4dforcing.test
  testoutput/4dhybrid.test
  testoutput/4dsaddlepoint.test
  testoutput/4dvar.alpha.test
  testoutput/4dvar.drgmresr.test
  testoutput/4dvar.dripcg.test
  testoutput/4dvar.dripcgqn.test
  testoutput/4dvar.drpcg.test
  testoutput/4dvar.drpcgqn.test
  testoutput/4dvar.drplanclmp.test
  testoutput/4dvar.drplanczos.test
  testoutput/4dvar.gmresr.test
  testoutput/4dvar.ipcg.test
  testoutput/4dvar.modbias.test
  testoutput/4dvar.obsbias.test
  testoutput/4dvar.pcg.test
  testoutput/4dvar.planczos.test
  testoutput/4dvar.rpcg.test
  testoutput/4dvar.rplanczos.test
  testoutput/forecast.test
  testoutput/genenspert.test
  testoutput/hofx.test
  testoutput/makeobs3d.test
  testoutput/makeobs4d.test
  testoutput/makeobs4d12h.test
  testoutput/makeobsbias.test
  testoutput/truth.test
)

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${l95_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for test input data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
foreach(FILENAME ${l95_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/testdata/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${l95_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Show in project
ecbuild_add_resources( TARGET   l95_test_scripts
                       SOURCES_PACK
                       ${l95_test_input}
                     )

#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

ecbuild_add_test( TARGET test_l95_truth
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_forecast.x testinput/truth.json"
                       testoutput/truth.test
                  DEPENDS l95_forecast.x )

ecbuild_add_test( TARGET test_l95_forecast
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_forecast.x testinput/forecast.json"
                        testoutput/forecast.test
                  DEPENDS l95_forecast.x )

ecbuild_add_test( TARGET test_l95_makeobs3d
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobs3d.json"
                       testoutput/makeobs3d.test
                  DEPENDS l95_makeobs.x )

ecbuild_add_test( TARGET test_l95_makeobs4d
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobs4d.json"
                       testoutput/makeobs4d.test
                  DEPENDS l95_makeobs.x )

ecbuild_add_test( TARGET test_l95_makeobsbias
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobsbias.json"
                       testoutput/makeobsbias.test
                  DEPENDS l95_makeobs.x )

ecbuild_add_test( TARGET test_l95_makeobs4d12h
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobs4d12h.json"
                       testoutput/makeobs4d12h.test
                  DEPENDS l95_makeobs.x )

# Test lorenz95 internal classes

ecbuild_add_test( TARGET test_l95_resolution
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/Resolution.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_stateL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/StateL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_tlmL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/TLML95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_ObsTable
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ObsTable.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsVec1D
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ObsVec1D.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_locsL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/LocsL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_gomL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/GomL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_observationL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ObservationL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsBias
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ObsBias.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_incrementL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/IncrementL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_localizationMatrixL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/LocalizationMatrixL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelBias
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ModelBias.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelBiasCorrection
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ModelBiasCorrection.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelBiasCovariance
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ModelBiasCovariance.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelL95
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ModelL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsBiasCorrection
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ObsBiasCorrection.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsBiasCovariance
                  BOOST
                  SOURCES lorenz95/test_main.cc lorenz95/ObsBiasCovariance.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.json"
                  LIBS    lorenz95 )

# Test interface classes with lorenz95

ecbuild_add_test( TARGET test_l95_geometry
                  BOOST
                  SOURCES executables/TestGeometry.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_locations
                  BOOST
                  SOURCES executables/TestLocations.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_variables
                  BOOST
                  SOURCES executables/TestVariables.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelaux
                  BOOST
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelauxincrement
                  BOOST
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelauxcovariance
                  BOOST
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_linearmodel
                  BOOST
                  SOURCES executables/TestLinearModel.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_observationspace
                  BOOST
                  SOURCES executables/TestObservationSpace.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsvector
                  BOOST
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obserrorcovar
                  BOOST
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsaux
                  BOOST
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsauxincrement
                  BOOST
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsauxcovariance
                  BOOST
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_localization
                  BOOST
                  SOURCES executables/TestLocalization.cc
                  ARGS "testinput/interfaces.json"
                  LIBS    lorenz95 )

# Test oops classes with lorenz95

#ecbuild_add_test( TARGET test_l95_op_obs.x
#                  BOOST
#                  SOURCES executables/TestSuiteOpObs.cc
#                  ARGS "testinput/test_op_obs.json"
#                  LIBS    lorenz95 )

# Test oops applications with lorenz95

ecbuild_add_test( TARGET test_l95_hofx
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_hofx.x testinput/hofx.json"
                       testoutput/hofx.test
                  DEPENDS l95_hofx.x )

ecbuild_add_test( TARGET test_l95_3dvar
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/3dvar.json"
                       testoutput/3dvar.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drgmresr
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drgmresr.json"
                       testoutput/4dvar.drgmresr.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_dripcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.dripcg.json"
                       testoutput/4dvar.dripcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_dripcgqn
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.dripcgqn.json"
                       testoutput/4dvar.dripcgqn.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drpcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drpcg.json"
                       testoutput/4dvar.drpcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drpcgqn
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drpcgqn.json"
                       testoutput/4dvar.drpcgqn.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drplanczos
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drplanczos.json"
                       testoutput/4dvar.drplanczos.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drplanclmp
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drplanclmp.json"
                       testoutput/4dvar.drplanclmp.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_gmresr
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.gmresr.json"
                       testoutput/4dvar.gmresr.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_ipcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.ipcg.json"
                       testoutput/4dvar.ipcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_pcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.pcg.json"
                       testoutput/4dvar.pcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_planczos
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.planczos.json"
                       testoutput/4dvar.planczos.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_rpcg
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.rpcg.json"
                       testoutput/4dvar.rpcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_rplanczos
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.rplanczos.json"
                       testoutput/4dvar.rplanczos.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_modbias
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.modbias.json"
                       testoutput/4dvar.modbias.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_forcing
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dforcing.json"
                       testoutput/4dforcing.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_saddlepoint
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dsaddlepoint.json"
                       testoutput/4dsaddlepoint.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_obsbias
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.obsbias.json"
                       testoutput/4dvar.obsbias.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_3dfgat
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/3dfgat.json"
                       testoutput/3dfgat.test
                  DEPENDS l95_4dvar.x )

if( ${OOPS_SENSITIVE_TESTS} )
ecbuild_add_test( TARGET test_l95_genenspert
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_genpert.x testinput/genenspert.json"
                       testoutput/genenspert.test
                  DEPENDS l95_genpert.x )

ecbuild_add_test( TARGET test_l95_alpha
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.alpha.json"
                       testoutput/4dvar.alpha.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_hybrid
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dhybrid.json"
                       testoutput/4dhybrid.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4densvar
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4densvar.json"
                       testoutput/4densvar.test
                  DEPENDS l95_4dvar.x )
endif()

