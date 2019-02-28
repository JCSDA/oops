list( APPEND l95_test_input
  testinput/3dvar.yaml
  testinput/3dfgat.yaml
  testinput/4dvar.alpha.yaml
  testinput/4dvar.drgmresr.yaml
  testinput/4dvar.dripcg.yaml
  testinput/4dvar.dripcgqn.yaml
  testinput/4dvar.drpcg.yaml
  testinput/4dvar.drpcgqn.yaml
  testinput/4dvar.drplanczos.yaml
  testinput/4dvar.drplanclmp.yaml
  testinput/4dvar.eda.pertobs.yaml
  testinput/4dvar.gmresr.yaml
  testinput/4dvar.ipcg.yaml
  testinput/4dvar.pcg.yaml
  testinput/4dvar.planczos.yaml
  testinput/4dvar.rpcg.yaml
  testinput/4dvar.rplanczos.yaml
  testinput/4dvar.modbias.yaml
  testinput/4dvar.obsbias.yaml
  testinput/4dforcing.yaml
  testinput/4dsaddlepoint.yaml
  testinput/4dhybrid.yaml
  testinput/4densvar.yaml
  testinput/enshofx.yaml
  testinput/forecast.yaml
  testinput/genenspert.yaml
  testinput/hofx.yaml
  testinput/interfaces.yaml
  testinput/makeobs3d.yaml
  testinput/makeobs4d.yaml
  testinput/makeobs4d12h.yaml
  testinput/makeobsbias.yaml
  testinput/test.yaml
  testinput/truth.yaml
  testinput/test_op_obs.yaml
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
  testoutput/enshofx.test
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
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_forecast.x testinput/truth.yaml"
                       testoutput/truth.test
                  DEPENDS l95_forecast.x )

ecbuild_add_test( TARGET test_l95_forecast
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_forecast.x testinput/forecast.yaml"
                        testoutput/forecast.test
                  DEPENDS l95_forecast.x )

ecbuild_add_test( TARGET test_l95_makeobs3d
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobs3d.yaml"
                       testoutput/makeobs3d.test
                  DEPENDS l95_makeobs.x )

ecbuild_add_test( TARGET test_l95_makeobs4d
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobs4d.yaml"
                       testoutput/makeobs4d.test
                  DEPENDS l95_makeobs.x )

ecbuild_add_test( TARGET test_l95_makeobsbias
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobsbias.yaml"
                       testoutput/makeobsbias.test
                  DEPENDS l95_makeobs.x )

ecbuild_add_test( TARGET test_l95_makeobs4d12h
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_makeobs.x testinput/makeobs4d12h.yaml"
                       testoutput/makeobs4d12h.test
                  DEPENDS l95_makeobs.x )

# Test lorenz95 internal classes

ecbuild_add_test( TARGET test_l95_resolution
                  SOURCES lorenz95/Resolution.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_stateL95
                  SOURCES lorenz95/StateL95.cc lorenz95/TestConfig.h
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  ARGS "testinput/test.yaml"
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_tlmL95
                  SOURCES lorenz95/TLML95.cc lorenz95/TestConfig.h
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  ARGS "testinput/test.yaml"
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_ObsTable
                  SOURCES lorenz95/ObsTable.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_obsVec1D
                  SOURCES lorenz95/ObsVec1D.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_locsL95
                  SOURCES lorenz95/LocsL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_gomL95
                  SOURCES lorenz95/GomL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_observationL95
                  SOURCES lorenz95/ObservationL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_obsBias
                  SOURCES lorenz95/ObsBias.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_incrementL95
                  SOURCES lorenz95/IncrementL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_localizationMatrixL95
                  SOURCES lorenz95/LocalizationMatrixL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_modelBias
                  SOURCES lorenz95/ModelBias.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_modelBiasCorrection
                  SOURCES lorenz95/ModelBiasCorrection.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_modelBiasCovariance
                  SOURCES lorenz95/ModelBiasCovariance.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_modelL95
                  SOURCES lorenz95/ModelL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_obsBiasCorrection
                  SOURCES lorenz95/ObsBiasCorrection.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 eckit)

ecbuild_add_test( TARGET test_l95_obsBiasCovariance
                  SOURCES lorenz95/ObsBiasCovariance.cc lorenz95/TestConfig.h
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  ARGS "testinput/test.yaml"
                  LIBS    lorenz95 eckit)

# Test interface classes with lorenz95

ecbuild_add_test( TARGET test_l95_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_state
                  SOURCES executables/TestState.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelaux
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_model
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_errorcovariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelauxincrement
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_modelauxcovariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_linearmodel
                  SOURCES executables/TestLinearModel.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_observationspace
                  SOURCES executables/TestObservationSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsvector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsoperator
                  SOURCES executables/TestObsOperator.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_linobsoperator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obserrorcovar
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsaux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsauxincrement
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_obsauxcovariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

ecbuild_add_test( TARGET test_l95_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    lorenz95 )

# Test oops applications with lorenz95

ecbuild_add_test( TARGET test_l95_hofx
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_hofx.x testinput/hofx.yaml"
                       testoutput/hofx.test
                  DEPENDS l95_hofx.x )

ecbuild_add_test( TARGET test_l95_3dvar
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/3dvar.yaml"
                       testoutput/3dvar.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drgmresr
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drgmresr.yaml"
                       testoutput/4dvar.drgmresr.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_dripcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.dripcg.yaml"
                       testoutput/4dvar.dripcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_dripcgqn
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.dripcgqn.yaml"
                       testoutput/4dvar.dripcgqn.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drpcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drpcg.yaml"
                       testoutput/4dvar.drpcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drpcgqn
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drpcgqn.yaml"
                       testoutput/4dvar.drpcgqn.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drplanczos
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drplanczos.yaml"
                       testoutput/4dvar.drplanczos.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_drplanclmp
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.drplanclmp.yaml"
                       testoutput/4dvar.drplanclmp.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_gmresr
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.gmresr.yaml"
                       testoutput/4dvar.gmresr.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_ipcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.ipcg.yaml"
                       testoutput/4dvar.ipcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_pcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.pcg.yaml"
                       testoutput/4dvar.pcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_planczos
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.planczos.yaml"
                       testoutput/4dvar.planczos.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_rpcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.rpcg.yaml"
                       testoutput/4dvar.rpcg.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_rplanczos
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.rplanczos.yaml"
                       testoutput/4dvar.rplanczos.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_modbias
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.modbias.yaml"
                       testoutput/4dvar.modbias.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_forcing
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dforcing.yaml"
                       testoutput/4dforcing.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_saddlepoint
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dsaddlepoint.yaml"
                       testoutput/4dsaddlepoint.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_obsbias
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.obsbias.yaml"
                       testoutput/4dvar.obsbias.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_3dfgat
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/3dfgat.yaml"
                       testoutput/3dfgat.test
                  DEPENDS l95_4dvar.x )

if( ${OOPS_SENSITIVE_TESTS} )
ecbuild_add_test( TARGET test_l95_genenspert
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_genpert.x testinput/genenspert.yaml"
                       testoutput/genenspert.test
                  DEPENDS l95_genpert.x )

ecbuild_add_test( TARGET test_l95_enshofx
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_enshofx.x testinput/enshofx.yaml"
                       testoutput/enshofx.test
                  DEPENDS l95_enshofx.x )

ecbuild_add_test( TARGET test_l95_alpha
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dvar.alpha.yaml"
                       testoutput/4dvar.alpha.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4dvar_hybrid
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4dhybrid.yaml"
                       testoutput/4dhybrid.test
                  DEPENDS l95_4dvar.x )

ecbuild_add_test( TARGET test_l95_4densvar
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/l95_4dvar.x testinput/4densvar.yaml"
                       testoutput/4densvar.test
                  DEPENDS l95_4dvar.x )
endif()

ecbuild_add_test( TARGET test_l95_4dvar_eda
                  MPI 4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/4dvar.eda.pertobs.yaml
                  DEPENDS l95_eda.x )

