include ( oops_add_test )

list( APPEND l95_test_input
  testinput/3dvar.yaml
  testinput/3dvar_qc.yaml
  testinput/3dvar_qc_iterations.yaml
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
  testinput/4dvar.allbiases.yaml
  testinput/4dforcing.yaml
  testinput/4dsaddlepoint.yaml
  testinput/4dhybrid.yaml
  testinput/4densvar.yaml
  testinput/diffstates.yaml
  testinput/eda.3dfgat.yaml
  testinput/eda.3dvar.yaml
  testinput/eda.4dvar.member3.yaml
  testinput/enshofx.yaml
  testinput/forecast.yaml
  testinput/forecast_pseudomodel.yaml
  testinput/forecast_identitymodel.yaml
  testinput/genenspert.yaml
  testinput/hofx.yaml
  testinput/hofx_nomodel.yaml
  testinput/identitymodel.yaml
  testinput/letkfOOPS.yaml
  testinput/letkf.yaml
  testinput/localhofx.yaml
  testinput/interfaces.yaml
  testinput/makeobs3d.yaml
  testinput/makeobs4d.yaml
  testinput/makeobs4d12h.yaml
  testinput/makeobsbias.yaml
  testinput/pseudomodel.yaml
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
  testoutput/3dvar_qc.test
  testoutput/3dvar_qc_iterations.test
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
  testoutput/4dvar.allbiases.test
  testoutput/4dvar.pcg.test
  testoutput/4dvar.planczos.test
  testoutput/4dvar.rpcg.test
  testoutput/4dvar.rplanczos.test
  testoutput/diffstates.test
  testoutput/enshofx.test
  testoutput/forecast.test
  testoutput/forecast_pseudomodel.test
  testoutput/forecast_identitymodel.test
  testoutput/genenspert.test
  testoutput/hofx.test
  testoutput/hofx.nomodel.test
  testoutput/letkf.test
  testoutput/letkf_oops.test
  testoutput/localhofx.test
  testoutput/makeobs3d.test
  testoutput/makeobs4d.test
  testoutput/makeobs4d12h.test
  testoutput/makeobsbias.test
  testoutput/truth.test
  testoutput/CTestCostData.txt.test
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
# internal tests
#####################################################################

ecbuild_add_test( TARGET test_l95_resolution
                  SOURCES lorenz95/Resolution.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_stateL95
                  SOURCES lorenz95/StateL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_tlmL95
                  SOURCES lorenz95/TLML95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_ObsTable
                  SOURCES lorenz95/ObsTable.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_obsVec1D
                  SOURCES lorenz95/ObsVec1D.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_locsL95
                  SOURCES lorenz95/LocsL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_gomL95
                  SOURCES lorenz95/GomL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_observationL95
                  SOURCES lorenz95/ObservationL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_incrementL95
                  SOURCES lorenz95/IncrementL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_localizationMatrixL95
                  SOURCES lorenz95/LocalizationMatrixL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_modelBias
                  SOURCES lorenz95/ModelBias.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_modelBiasCorrection
                  SOURCES lorenz95/ModelBiasCorrection.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_modelBiasCovariance
                  SOURCES lorenz95/ModelBiasCovariance.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_modelL95
                  SOURCES lorenz95/ModelL95.cc lorenz95/TestConfig.h
                  ARGS "testinput/test.yaml"
                  LIBS lorenz95 )


#####################################################################
# truth and makeobs4d are required for interface tests
#####################################################################

oops_add_test( TESTNAME truth
               MODELNAME l95
               YAMLNAME testinput/truth.yaml
               EXENAME  l95_forecast.x )

oops_add_test( TESTNAME makeobs4d
               MODELNAME l95
               YAMLNAME testinput/makeobs4d.yaml
               EXENAME l95_makeobs.x
               TEST_DEPENDS test_l95_truth)


#####################################################################
# interface tests
#####################################################################

ecbuild_add_test( TARGET test_l95_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_geometry_iterator
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_state
                  SOURCES executables/TestState.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_modelaux
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_model
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_pseudomodel
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/pseudomodel.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_identitymodel
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/identitymodel.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_errorcovariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_modelauxincrement
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_modelauxcovariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_linearmodel
                  SOURCES executables/TestLinearModel.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsspace
                  SOURCES executables/TestObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_localobsspace
                  SOURCES executables/TestLocalObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsvector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_departures_ensemble
                  SOURCES executables/TestDeparturesEnsemble.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_departures
                  SOURCES executables/TestDepartures.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsoperator
                  SOURCES executables/TestObsOperator.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_linobsoperator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obserrorcovar
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsaux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsauxcovariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsauxincrement
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )


#####################################################################
# forecast-related tests
#####################################################################

oops_add_test( TESTNAME forecast
               MODELNAME l95
               YAMLNAME testinput/forecast.yaml
               EXENAME l95_forecast.x )

oops_add_test( TESTNAME forecast_pseudomodel
               MODELNAME l95
               YAMLNAME testinput/forecast_pseudomodel.yaml
               EXENAME l95_forecast.x )

oops_add_test( TESTNAME forecast_identitymodel
               MODELNAME l95
               YAMLNAME testinput/forecast_identitymodel.yaml
               EXENAME l95_forecast.x )


#####################################################################
# obs-related tests
#####################################################################

oops_add_test( TESTNAME makeobs3d
               MODELNAME l95
               YAMLNAME testinput/makeobs3d.yaml
               EXENAME l95_makeobs.x
               TEST_DEPENDS test_l95_truth ) 

oops_add_test( TESTNAME makeobsbias
               MODELNAME l95
               YAMLNAME testinput/makeobsbias.yaml
               EXENAME l95_makeobs.x
               TEST_DEPENDS test_l95_truth ) 

oops_add_test( TESTNAME makeobs4d12h
               MODELNAME l95
               YAMLNAME testinput/makeobs4d12h.yaml
               EXENAME l95_makeobs.x
               TEST_DEPENDS test_l95_truth ) 

oops_add_test( TESTNAME hofx
               MODELNAME l95
               YAMLNAME testinput/hofx.yaml
               EXENAME l95_hofx.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME hofx.nomodel
               MODELNAME l95
               YAMLNAME testinput/hofx_nomodel.yaml
               EXENAME l95_hofx_nomodel.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME localhofx
               MODELNAME l95
               YAMLNAME testinput/localhofx.yaml
               EXENAME l95_localhofx.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d ) 


#####################################################################
# ensemble-related tests
#####################################################################

oops_add_test( TESTNAME genenspert
               MODELNAME l95
               YAMLNAME testinput/genenspert.yaml
               EXENAME l95_genpert.x )

ecbuild_add_test( TARGET test_l95_enshofx
                  MPI 10
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_enshofx.x
                  ARGS testinput/enshofx.yaml testoutput/enshofx.log.out
                  DEPENDS l95_enshofx.x
                  TEST_DEPENDS test_l95_genenspert test_l95_makeobs4d )


#####################################################################
# 3d variational tests
#####################################################################

oops_add_test( TESTNAME 3dvar
               MODELNAME l95
               YAMLNAME testinput/3dvar.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob3d ) 

oops_add_test( TESTNAME 3dvar_qc
               MODELNAME l95
               YAMLNAME testinput/3dvar_qc.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob3d )

oops_add_test( TESTNAME 3dvar_qc_iterations
               MODELNAME l95
               YAMLNAME testinput/3dvar_qc_iterations.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob3d )

oops_add_test( TESTNAME 3dfgat
               MODELNAME l95
               YAMLNAME testinput/3dfgat.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob3d ) 


#####################################################################
# 4d variational tests
#####################################################################

oops_add_test( TESTNAME 4dvar.drgmresr
               MODELNAME l95
               YAMLNAME testinput/4dvar.drgmresr.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d ) 

oops_add_test( TESTNAME 4dvar.dripcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.dripcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d ) 

oops_add_test( TESTNAME 4dvar.dripcgqn
               MODELNAME l95
               YAMLNAME testinput/4dvar.dripcgqn.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d ) 

oops_add_test( TESTNAME 4dvar.drpcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.drpcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.drpcgqn
               MODELNAME l95
               YAMLNAME testinput/4dvar.drpcgqn.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.drplanczos
               MODELNAME l95
               YAMLNAME testinput/4dvar.drplanczos.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.drplanclmp
               MODELNAME l95
               YAMLNAME testinput/4dvar.drplanclmp.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.gmresr
               MODELNAME l95
               YAMLNAME testinput/4dvar.gmresr.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.ipcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.ipcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.pcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.pcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.planczos
               MODELNAME l95
               YAMLNAME testinput/4dvar.planczos.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.rpcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.rpcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.rplanczos
               MODELNAME l95
               YAMLNAME testinput/4dvar.rplanczos.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.modbias
               MODELNAME l95
               YAMLNAME testinput/4dvar.modbias.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dforcing
               MODELNAME l95
               YAMLNAME testinput/4dforcing.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dsaddlepoint
               MODELNAME l95
               YAMLNAME testinput/4dsaddlepoint.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d )

oops_add_test( TESTNAME 4dvar.obsbias
               MODELNAME l95
               YAMLNAME testinput/4dvar.obsbias.yaml
               EXENAME l95_4dvar.x
               COMPARE
               TEST_DEPENDS test_l95_forecast test_l95_makeobsbias )
#--------------------------------------------------------------------

oops_add_test( TESTNAME 4dvar.allbiases
               MODELNAME l95
               YAMLNAME testinput/4dvar.allbiases.yaml
               EXENAME l95_4dvar.x
               COMPARE
               TEST_DEPENDS test_l95_forecast test_l95_makeobsbias )
#--------------------------------------------------------------------

oops_add_test( TESTNAME 4dvar.alpha
               MODELNAME l95
               YAMLNAME testinput/4dvar.alpha.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d test_l95_genenspert)

oops_add_test( TESTNAME 4dhybrid
               MODELNAME l95
               YAMLNAME testinput/4dhybrid.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d test_l95_genenspert )

oops_add_test( TESTNAME 4densvar
               MODELNAME l95
               YAMLNAME testinput/4densvar.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeob4d12h test_l95_genenspert )


#####################################################################
# EDA tests
#####################################################################

ecbuild_add_test( TARGET test_l95_eda_3dfgat
                  MPI 4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/eda.3dfgat.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_eda_3dvar
                  MPI 4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/eda.3dvar.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_eda_4dvar_member3
                  MPI 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/eda.4dvar.member3.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )


#####################################################################
# state-related tests
#####################################################################

oops_add_test( TESTNAME diffstates
               MODELNAME l95
               YAMLNAME testinput/diffstates.yaml
               EXENAME l95_diffstates.x
               TEST_DEPENDS test_l95_eda_3dvar  test_l95_eda_4dvar_member3 )


#####################################################################
# LETKF tests
#####################################################################

oops_add_test( TESTNAME letkf_oops
               MODELNAME l95
               YAMLNAME testinput/letkfOOPS.yaml
               EXENAME l95_letkf.x
               COMPARE
               TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

oops_add_test( TESTNAME letkf
               MODELNAME l95
               YAMLNAME testinput/letkf.yaml
               EXENAME l95_letkf.x
               TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )
