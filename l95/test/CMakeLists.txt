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
  testinput/eda.3dfgat.yaml
  testinput/eda.3dvar.yaml
  testinput/eda.4dvar.yaml
  testinput/enshofx.yaml
  testinput/forecast.yaml
  testinput/genenspert.yaml
  testinput/hofx.yaml
  testinput/letkf.yaml
  testinput/localhofx.yaml
  testinput/interfaces.yaml
  testinput/makeobs3d.yaml
  testinput/makeobs4d.yaml
  testinput/makeobs4d12h.yaml
  testinput/makeobsbias.yaml
  testinput/test.yaml
  testinput/truth.yaml
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
  testoutput/letkf.test
  testoutput/localhofx.test
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
# internal tests
#####################################################################

ecbuild_add_test( TARGET test_l95_resolution
                  SOURCES lorenz95/Resolution.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_stateL95
                  SOURCES lorenz95/StateL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_tlmL95
                  SOURCES lorenz95/TLML95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_ObsTable
                  SOURCES lorenz95/ObsTable.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obsVec1D
                  SOURCES lorenz95/ObsVec1D.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_locsL95
                  SOURCES lorenz95/LocsL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_gomL95
                  SOURCES lorenz95/GomL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_observationL95
                  SOURCES lorenz95/ObservationL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_incrementL95
                  SOURCES lorenz95/IncrementL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_localizationMatrixL95
                  SOURCES lorenz95/LocalizationMatrixL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_modelBias
                  SOURCES lorenz95/ModelBias.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_modelBiasCorrection
                  SOURCES lorenz95/ModelBiasCorrection.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_modelBiasCovariance
                  SOURCES lorenz95/ModelBiasCovariance.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_modelL95
                  SOURCES lorenz95/ModelL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 eckit )

#####################################################################
# truth and makeobs4d are required for interface tests
#####################################################################

ecbuild_add_test( TARGET test_l95_truth_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_forecast.x
                  ARGS testinput/truth.yaml testoutput/truth.log.out
                  DEPENDS l95_forecast.x )

ecbuild_add_test( TARGET test_l95_truth_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/truth.log.out testoutput/truth.test
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_makeobs4d_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_makeobs.x
                  ARGS testinput/makeobs4d.yaml testoutput/makeobs4d.log.out
                  DEPENDS l95_makeobs.x
                  TEST_DEPENDS test_l95_truth_run )

ecbuild_add_test( TARGET test_l95_makeobs4d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/makeobs4d.log.out testoutput/makeobs4d.test
                  TEST_DEPENDS test_l95_makeobs4d_run )

#####################################################################
# interface tests
#####################################################################


ecbuild_add_test( TARGET test_l95_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_geometry_iterator
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_state
                  SOURCES executables/TestState.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_modelaux
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_model
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_errorcovariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_modelauxincrement
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_modelauxcovariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_linearmodel
                  SOURCES executables/TestLinearModel.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obsspace
                  SOURCES executables/TestObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_localobsspace
                  SOURCES executables/TestLocalObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obsvector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obsoperator
                  SOURCES executables/TestObsOperator.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_linobsoperator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obserrorcovar
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obsaux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obsauxcovariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_obsauxincrement
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth_run )

#####################################################################
# forecast-related tests
#####################################################################

ecbuild_add_test( TARGET test_l95_forecast_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_forecast.x
                  ARGS testinput/forecast.yaml testoutput/forecast.log.out
                  DEPENDS l95_forecast.x )

ecbuild_add_test( TARGET test_l95_forecast_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/forecast.log.out testoutput/forecast.test
                  TEST_DEPENDS test_l95_forecast_run )

#####################################################################
# obs-related tests
#####################################################################

ecbuild_add_test( TARGET test_l95_makeobs3d_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_makeobs.x
                  ARGS testinput/makeobs3d.yaml testoutput/makeobs3d.log.out
                  DEPENDS l95_makeobs.x
                  TEST_DEPENDS test_l95_truth_run )

ecbuild_add_test( TARGET test_l95_makeobs3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/makeobs3d.log.out testoutput/makeobs3d.test
                  TEST_DEPENDS test_l95_makeobs3d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_makeobsbias_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_makeobs.x
                  ARGS testinput/makeobsbias.yaml testoutput/makeobsbias.log.out
                  DEPENDS l95_makeobs.x
                  TEST_DEPENDS test_l95_truth_run )

ecbuild_add_test( TARGET test_l95_makeobsbias_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/makeobsbias.log.out testoutput/makeobsbias.test
                  TEST_DEPENDS test_l95_makeobsbias_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_makeobs4d12h_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_makeobs.x
                  ARGS testinput/makeobs4d12h.yaml testoutput/makeobs4d12h.log.out
                  DEPENDS l95_makeobs.x
                  TEST_DEPENDS test_l95_truth_run )

ecbuild_add_test( TARGET test_l95_makeobs4d12h_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/makeobs4d12h.log.out testoutput/makeobs4d12h.test
                  TEST_DEPENDS test_l95_makeobs4d12h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_hofx_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_hofx.x
                  ARGS testinput/hofx.yaml testoutput/hofx.log.out
                  DEPENDS l95_hofx.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_hofx_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/hofx.log.out testoutput/hofx.test
                  TEST_DEPENDS test_l95_hofx_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_localhofx_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_localhofx.x
                  ARGS testinput/localhofx.yaml testoutput/localhofx.log.out
                  DEPENDS l95_localhofx.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_localhofx_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/localhofx.log.out testoutput/localhofx.test
                  TEST_DEPENDS test_l95_localhofx_run )

#####################################################################
# ensemble-related tests
#####################################################################

ecbuild_add_test( TARGET test_l95_genenspert_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_genpert.x
                  ARGS testinput/genenspert.yaml testoutput/genenspert.log.out
                  DEPENDS l95_genpert.x )

ecbuild_add_test( TARGET test_l95_genenspert_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/genenspert.log.out testoutput/genenspert.test
                  TEST_DEPENDS test_l95_genenspert_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_enshofx_run
                  MPI 10
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_enshofx.x
                  ARGS testinput/enshofx.yaml testoutput/enshofx.log.out
                  DEPENDS l95_enshofx.x
                  TEST_DEPENDS test_l95_genenspert_run test_l95_makeobs4d_run )

#ecbuild_add_test( TARGET test_l95_enshofx_compare
#                  TYPE SCRIPT
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
#                  ARGS testoutput/enshofx.log.out testoutput/enshofx.test
#                  TEST_DEPENDS test_l95_enshofx_run )

#####################################################################
# 3d variational tests
#####################################################################

ecbuild_add_test( TARGET test_l95_3dvar_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/3dvar.yaml testoutput/3dvar.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs3d_run )

ecbuild_add_test( TARGET test_l95_3dvar_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar.log.out testoutput/3dvar.test
                  TEST_DEPENDS test_l95_3dvar_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_3dfgat_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/3dfgat.yaml testoutput/3dfgat.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs3d_run )

ecbuild_add_test( TARGET test_l95_3dfgat_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dfgat.log.out testoutput/3dfgat.test
                  TEST_DEPENDS test_l95_3dfgat_run )

#####################################################################
# 4d variational tests
#####################################################################

ecbuild_add_test( TARGET test_l95_4dvar_drgmresr_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.drgmresr.yaml testoutput/4dvar.drgmresr.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_drgmresr_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.drgmresr.log.out testoutput/4dvar.drgmresr.test
                  TEST_DEPENDS test_l95_4dvar_drgmresr_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_dripcg_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.dripcg.yaml testoutput/4dvar.dripcg.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_dripcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.dripcg.log.out testoutput/4dvar.dripcg.test
                  TEST_DEPENDS test_l95_4dvar_dripcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_dripcgqn_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.dripcgqn.yaml testoutput/4dvar.dripcgqn.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_dripcgqn_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.dripcgqn.log.out testoutput/4dvar.dripcgqn.test
                  TEST_DEPENDS test_l95_4dvar_dripcgqn_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_drpcg_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.drpcg.yaml testoutput/4dvar.drpcg.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_drpcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.drpcg.log.out testoutput/4dvar.drpcg.test
                  TEST_DEPENDS test_l95_4dvar_drpcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_drpcgqn_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.drpcgqn.yaml testoutput/4dvar.drpcgqn.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_drpcgqn_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.drpcgqn.log.out testoutput/4dvar.drpcgqn.test
                  TEST_DEPENDS test_l95_4dvar_drpcgqn_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_drplanczos_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.drplanczos.yaml testoutput/4dvar.drplanczos.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_drplanczos_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.drplanczos.log.out testoutput/4dvar.drplanczos.test
                  TEST_DEPENDS test_l95_4dvar_drplanczos_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_drplanclmp_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.drplanclmp.yaml testoutput/4dvar.drplanclmp.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_drplanclmp_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.drplanclmp.log.out testoutput/4dvar.drplanclmp.test
                  TEST_DEPENDS test_l95_4dvar_drplanclmp_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_gmresr_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.gmresr.yaml testoutput/4dvar.gmresr.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_gmresr_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.gmresr.log.out testoutput/4dvar.gmresr.test
                  TEST_DEPENDS test_l95_4dvar_gmresr_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_ipcg_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.ipcg.yaml testoutput/4dvar.ipcg.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_ipcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.ipcg.log.out testoutput/4dvar.ipcg.test
                  TEST_DEPENDS test_l95_4dvar_ipcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_pcg_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.pcg.yaml testoutput/4dvar.pcg.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_pcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.pcg.log.out testoutput/4dvar.pcg.test
                  TEST_DEPENDS test_l95_4dvar_pcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_planczos_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.planczos.yaml testoutput/4dvar.planczos.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_planczos_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.planczos.log.out testoutput/4dvar.planczos.test
                  TEST_DEPENDS test_l95_4dvar_planczos_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_rpcg_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.rpcg.yaml testoutput/4dvar.rpcg.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_rpcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.rpcg.log.out testoutput/4dvar.rpcg.test
                  TEST_DEPENDS test_l95_4dvar_rpcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_rplanczos_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.rplanczos.yaml testoutput/4dvar.rplanczos.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_rplanczos_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.rplanczos.log.out testoutput/4dvar.rplanczos.test
                  TEST_DEPENDS test_l95_4dvar_rplanczos_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_modbias_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.modbias.yaml testoutput/4dvar.modbias.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_modbias_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.modbias.log.out testoutput/4dvar.modbias.test
                  TEST_DEPENDS test_l95_4dvar_modbias_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_forcing_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dforcing.yaml testoutput/4dforcing.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_forcing_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dforcing.log.out testoutput/4dforcing.test
                  TEST_DEPENDS test_l95_4dvar_forcing_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_saddlepoint_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dsaddlepoint.yaml testoutput/4dsaddlepoint.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_4dvar_saddlepoint_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dsaddlepoint.log.out testoutput/4dsaddlepoint.test
                  TEST_DEPENDS test_l95_4dvar_saddlepoint_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_obsbias_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.obsbias.yaml testoutput/4dvar.obsbias.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobsbias_run )

ecbuild_add_test( TARGET test_l95_4dvar_obsbias_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.obsbias.log.out testoutput/4dvar.obsbias.test
                  TEST_DEPENDS test_l95_4dvar_obsbias_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_alpha_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dvar.alpha.yaml testoutput/4dvar.alpha.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run test_l95_genenspert_run)

ecbuild_add_test( TARGET test_l95_4dvar_alpha_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar.alpha.log.out testoutput/4dvar.alpha.test
                  TEST_DEPENDS test_l95_4dvar_alpha_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4dvar_hybrid_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4dhybrid.yaml testoutput/4dhybrid.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d_run test_l95_genenspert_run)

ecbuild_add_test( TARGET test_l95_4dvar_hybrid_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dhybrid.log.out testoutput/4dhybrid.test
                  TEST_DEPENDS test_l95_4dvar_hybrid_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_l95_4densvar_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_4dvar.x
                  ARGS testinput/4densvar.yaml testoutput/4densvar.log.out
                  DEPENDS l95_4dvar.x
                  TEST_DEPENDS test_l95_forecast_run test_l95_makeobs4d12h_run test_l95_genenspert_run)

ecbuild_add_test( TARGET test_l95_4densvar_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4densvar.log.out testoutput/4densvar.test
                  TEST_DEPENDS test_l95_4densvar_run )

#####################################################################
# EDA tests
#####################################################################

ecbuild_add_test( TARGET test_l95_eda_3dfgat
                  MPI 4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/eda.3dfgat.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast_run  test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_eda_3dvar
                  MPI 4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/eda.3dvar.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast_run  test_l95_makeobs4d_run )

ecbuild_add_test( TARGET test_l95_eda_4dvar
                  MPI 4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/eda.4dvar.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast_run  test_l95_makeobs4d_run )

#####################################################################
# LETKF tests
#####################################################################

ecbuild_add_test( TARGET test_l95_letkf_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_letkf.x
                  ARGS testinput/letkf.yaml testoutput/letkf.log.out
                  DEPENDS l95_letkf.x
                  TEST_DEPENDS test_l95_makeobs3d_run test_l95_genenspert_run)

ecbuild_add_test( TARGET test_l95_letkf_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/letkf.log.out testoutput/letkf.test
                  TEST_DEPENDS test_l95_letkf_run )
