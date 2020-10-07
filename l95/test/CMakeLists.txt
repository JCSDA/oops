include ( oops_add_test )

list( APPEND l95_test_input
  testinput/3dvar.yaml
  testinput/3dvar_qc.yaml
  testinput/3dvar_qc_obserr.yaml
  testinput/3dvar_qc_iterations.yaml
  testinput/3dfgat.yaml
  testinput/4densvar.yaml
  testinput/4densvar.hybrid.yaml
  testinput/4dforcing.yaml
  testinput/4dsaddlepoint.yaml
  testinput/4dvar.alpha.yaml
  testinput/4dvar.drgmresr.yaml
  testinput/4dvar.dripcg.yaml
  testinput/4dvar.dripcgqn.yaml
  testinput/4dvar.drpcg.yaml
  testinput/4dvar.drpcgqn.yaml
  testinput/4dvar.drplanczos.yaml
  testinput/4dvar.drplanclmp.yaml
  testinput/4dvar.fgmres.yaml
  testinput/4dvar.gmresr.yaml
  testinput/4dvar.hybrid.yaml
  testinput/4dvar.ipcg.yaml
  testinput/4dvar.lbgmresr.yaml
  testinput/4dvar.pcg.yaml
  testinput/4dvar.planczos.yaml
  testinput/4dvar.rpcg.yaml
  testinput/4dvar.rplanczos.yaml
  testinput/4dvar.minres.yaml
  testinput/4dvar.modbias.yaml
  testinput/4dvar.obsbias.yaml
  testinput/4dvar.allbiases.yaml
  testinput/addincrement.yaml
  testinput/addincrement_scaled.yaml
  testinput/diffstates.yaml
  testinput/eda.3dfgat.1.yaml
  testinput/eda.3dfgat.2.yaml
  testinput/eda.3dfgat.3.yaml
  testinput/eda.3dfgat.4.yaml
  testinput/eda.3dfgat.yaml
  testinput/eda.3dvar.1.yaml
  testinput/eda.3dvar.2.yaml
  testinput/eda.3dvar.3.yaml
  testinput/eda.3dvar.4.yaml
  testinput/eda.3dvar.yaml
  testinput/eda.4dvar.3.yaml
  testinput/eda.4dvar.yaml
  testinput/enshofx.1.yaml
  testinput/enshofx.2.yaml
  testinput/enshofx.3.yaml
  testinput/enshofx.4.yaml
  testinput/enshofx.yaml
  testinput/forecast.yaml
  testinput/forecast_pseudomodel.yaml
  testinput/forecast_identitymodel.yaml
  testinput/genenspert.yaml
  testinput/getkf.yaml
  testinput/hofx.yaml
  testinput/hofx_nomodel.yaml
  testinput/identitymodel.yaml
  testinput/letkfGSI.yaml
  testinput/letkf.yaml
  testinput/interfaces.yaml
  testinput/linearmodelfactory.yaml
  testinput/makeobs3d.yaml
  testinput/makeobs4d.yaml
  testinput/makeobs4d12h.yaml
  testinput/makeobsbias.yaml
  testinput/makeobspert.yaml
  testinput/pseudomodel.yaml
  testinput/truth.yaml
  testinput/simplifiedl95_DRGMRESR.yaml
  testinput/simplifiedl95_DRIPCG.yaml
  testinput/simplifiedl95_DRPCG.yaml
  testinput/simplifiedl95_DRPFOM.yaml
  testinput/simplifiedl95_DRPLanczos.yaml
  testinput/simplifiedl95_FGMRES.yaml
  testinput/simplifiedl95_GMRESR.yaml
  testinput/simplifiedl95_IPCG.yaml
  testinput/simplifiedl95_LBGMRESR.yaml
  testinput/simplifiedl95_MINRES.yaml
  testinput/simplifiedl95_PCG.yaml
  testinput/simplifiedl95_PLanczos.yaml
  testinput/simplifiedl95_RPCG.yaml
  testinput/simplifiedl95_RPLanczos.yaml
)

list( APPEND l95_test_data
  test.an.2010-01-01T00:00:00Z.l95
  truth.an.2010-01-01T00:00:00Z.l95
  truth.fc.2010-01-01T00:00:00Z.PT3H
  truth.2010-01-02T00:00:00Z.obt
  simplifiedl95.fc
  simplifiedl95.truth.obt
)

list( APPEND l95_testoutput
  testoutput/3dvar.test
  testoutput/3dvar_qc.test
  testoutput/3dvar_qc_obserr.test
  testoutput/3dvar_qc_iterations.test
  testoutput/3dfgat.test
  testoutput/4densvar.test
  testoutput/4densvar.hybrid.test
  testoutput/4dforcing.test
  testoutput/4dsaddlepoint.test
  testoutput/4dvar.alpha.test
  testoutput/4dvar.drgmresr.test
  testoutput/4dvar.dripcg.test
  testoutput/4dvar.dripcgqn.test
  testoutput/4dvar.drpcg.test
  testoutput/4dvar.drpcgqn.test
  testoutput/4dvar.drplanclmp.test
  testoutput/4dvar.drplanczos.test
  testoutput/4dvar.fgmres.test
  testoutput/4dvar.gmresr.test
  testoutput/4dvar.hybrid.test
  testoutput/4dvar.ipcg.test
  testoutput/4dvar.lbgmresr.test
  testoutput/4dvar.modbias.test
  testoutput/4dvar.obsbias.test
  testoutput/4dvar.allbiases.test
  testoutput/4dvar.pcg.test
  testoutput/4dvar.planczos.test
  testoutput/4dvar.rpcg.test
  testoutput/4dvar.rplanczos.test
  testoutput/4dvar.minres.test
  testoutput/diffstates.test
  testoutput/addincrement.test
  testoutput/addincrement_scaled.test
  testoutput/enshofx.test
  testoutput/forecast.test
  testoutput/forecast_pseudomodel.test
  testoutput/forecast_identitymodel.test
  testoutput/genenspert.test
  testoutput/getkf.test
  testoutput/hofx.test
  testoutput/hofx.nomodel.test
  testoutput/letkf.test
  testoutput/letkf_gsi.test
  testoutput/localhofx.test
  testoutput/makeobs3d.test
  testoutput/makeobs4d.test
  testoutput/makeobs4d12h.test
  testoutput/makeobsbias.test
  testoutput/makeobspert.test
  testoutput/truth.test
  testoutput/CTestCostData.txt.test
  testoutput/simplifiedl95_DRGMRESR.test
  testoutput/simplifiedl95_DRIPCG.test
  testoutput/simplifiedl95_DRPCG.test
  testoutput/simplifiedl95_DRPFOM.test
  testoutput/simplifiedl95_DRPLanczos.test
  testoutput/simplifiedl95_FGMRES.test
  testoutput/simplifiedl95_GMRESR.test
  testoutput/simplifiedl95_IPCG.test
  testoutput/simplifiedl95_LBGMRESR.test
  testoutput/simplifiedl95_MINRES.test
  testoutput/simplifiedl95_PCG.test
  testoutput/simplifiedl95_PLanczos.test
  testoutput/simplifiedl95_RPCG.test
  testoutput/simplifiedl95_RPLanczos.test
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
# truth and makeobs4d are required for interface tests
#####################################################################

oops_add_test( TESTNAME truth
               MODELNAME l95
               YAMLNAME testinput/truth.yaml
               EXENAME  l95_forecast.x )

oops_add_test( TESTNAME makeobs4d
               MODELNAME l95
               YAMLNAME testinput/makeobs4d.yaml
               EXENAME l95_hofx.x
               TEST_DEPENDS test_l95_truth )


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

ecbuild_add_test( TARGET test_l95_getvalues
                  SOURCES executables/TestGetValues.cc
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

ecbuild_add_test( TARGET test_l95_lineargetvalues
                  SOURCES executables/TestLinearGetValues.cc
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
               EXENAME l95_hofx.x
               TEST_DEPENDS test_l95_truth )

oops_add_test( TESTNAME makeobsbias
               MODELNAME l95
               YAMLNAME testinput/makeobsbias.yaml
               EXENAME l95_hofx.x
               TEST_DEPENDS test_l95_truth )

oops_add_test( TESTNAME makeobs4d12h
               MODELNAME l95
               YAMLNAME testinput/makeobs4d12h.yaml
               EXENAME l95_hofx.x
               TEST_DEPENDS test_l95_truth )

oops_add_test( TESTNAME makeobspert
               MODELNAME l95
               YAMLNAME testinput/makeobspert.yaml
               EXENAME l95_hofx.x
               TEST_DEPENDS test_l95_truth )

oops_add_test( TESTNAME hofx
               MODELNAME l95
               YAMLNAME testinput/hofx.yaml
               EXENAME l95_hofx.x
               TEST_DEPENDS test_l95_truth test_l95_makeobs4d )

oops_add_test( TESTNAME hofx.nomodel
               MODELNAME l95
               YAMLNAME testinput/hofx_nomodel.yaml
               EXENAME l95_hofx_nomodel.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

#####################################################################
# ensemble-related tests
#####################################################################

oops_add_test( TESTNAME genenspert
               MODELNAME l95
               YAMLNAME testinput/genenspert.yaml
               EXENAME l95_genpert.x )

ecbuild_add_test( TARGET test_l95_enshofx
                  MPI 4
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
               TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

oops_add_test( TESTNAME 3dvar_qc
               MODELNAME l95
               YAMLNAME testinput/3dvar_qc.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

oops_add_test( TESTNAME 3dvar_qc_obserr
               MODELNAME l95
               YAMLNAME testinput/3dvar_qc_obserr.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

oops_add_test( TESTNAME 3dvar_qc_iterations
               MODELNAME l95
               YAMLNAME testinput/3dvar_qc_iterations.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

oops_add_test( TESTNAME 3dfgat
               MODELNAME l95
               YAMLNAME testinput/3dfgat.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )


#####################################################################
# 4d variational tests
#####################################################################

oops_add_test( TESTNAME 4dvar.drgmresr
               MODELNAME l95
               YAMLNAME testinput/4dvar.drgmresr.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.dripcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.dripcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.dripcgqn
               MODELNAME l95
               YAMLNAME testinput/4dvar.dripcgqn.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.drpcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.drpcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.drpcgqn
               MODELNAME l95
               YAMLNAME testinput/4dvar.drpcgqn.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.drplanczos
               MODELNAME l95
               YAMLNAME testinput/4dvar.drplanczos.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.drplanclmp
               MODELNAME l95
               YAMLNAME testinput/4dvar.drplanclmp.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.fgmres
               MODELNAME l95
               YAMLNAME testinput/4dvar.fgmres.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.gmresr
               MODELNAME l95
               YAMLNAME testinput/4dvar.gmresr.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.ipcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.ipcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.lbgmresr
               MODELNAME l95
               YAMLNAME testinput/4dvar.lbgmresr.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.pcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.pcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.planczos
               MODELNAME l95
               YAMLNAME testinput/4dvar.planczos.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.rpcg
               MODELNAME l95
               YAMLNAME testinput/4dvar.rpcg.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.rplanczos
               MODELNAME l95
               YAMLNAME testinput/4dvar.rplanczos.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.minres
               MODELNAME l95
               YAMLNAME testinput/4dvar.minres.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dvar.modbias
               MODELNAME l95
               YAMLNAME testinput/4dvar.modbias.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

#oops_add_test( TESTNAME 4dforcing
#               MODELNAME l95
#               YAMLNAME testinput/4dforcing.yaml
#               EXENAME l95_4dvar.x
#               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

oops_add_test( TESTNAME 4dsaddlepoint
               MPI 2
               MODELNAME l95
               YAMLNAME testinput/4dsaddlepoint.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

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
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d test_l95_genenspert )

oops_add_test( TESTNAME 4dvar.hybrid
               MODELNAME l95
               YAMLNAME testinput/4dvar.hybrid.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d test_l95_genenspert )

oops_add_test( TESTNAME 4densvar
               MPI 9
               MODELNAME l95
               YAMLNAME testinput/4densvar.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d12h test_l95_genenspert )

oops_add_test( TESTNAME 4densvar.hybrid
               MPI 9
               MODELNAME l95
               YAMLNAME testinput/4densvar.hybrid.yaml
               EXENAME l95_4dvar.x
               TEST_DEPENDS test_l95_forecast test_l95_makeobs4d test_l95_genenspert )


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

ecbuild_add_test( TARGET test_l95_eda_4dvar
                  MPI 1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/l95_eda.x
                  ARGS testinput/eda.4dvar.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )


#####################################################################
# state-related tests
#####################################################################

oops_add_test( TESTNAME diffstates
               MODELNAME l95
               YAMLNAME testinput/diffstates.yaml
               EXENAME l95_diffstates.x
               TEST_DEPENDS test_l95_eda_3dvar test_l95_eda_4dvar )

oops_add_test( TESTNAME addincrement
               MODELNAME l95
               YAMLNAME testinput/addincrement.yaml
               EXENAME l95_addincrement.x
               TEST_DEPENDS test_l95_diffstates )

oops_add_test( TESTNAME addincrement_scaled
               MODELNAME l95
               YAMLNAME testinput/addincrement_scaled.yaml
               EXENAME l95_addincrement.x
               TEST_DEPENDS test_l95_diffstates )


#####################################################################
# LETKF tests
#####################################################################

oops_add_test( TESTNAME letkf
               MODELNAME l95
               YAMLNAME testinput/letkf.yaml
               EXENAME l95_letkf.x
               TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

oops_add_test( TESTNAME letkf_gsi
               MODELNAME l95
               YAMLNAME testinput/letkfGSI.yaml
               EXENAME l95_letkf.x
               CTOL 1e-4
               TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

oops_add_test( TESTNAME getkf
               MODELNAME l95
               YAMLNAME testinput/getkf.yaml
               EXENAME l95_letkf.x
               CTOL 1e-4
               TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

#####################################################################
# tests of factory and parameters classes
#####################################################################

ecbuild_add_test( TARGET test_l95_linearmodelfactory
                  SOURCES lorenz95/LinearModelFactory.cc
                  ARGS "testinput/linearmodelfactory.yaml"
                  LIBS lorenz95 )

#####################################################################
# tests with simplified L95
#####################################################################

oops_add_test( TESTNAME simplifiedl95_DRGMRESR
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_DRGMRESR.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_DRIPCG
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_DRIPCG.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_DRPCG
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_DRPCG.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_DRPFOM
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_DRPFOM.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_DRPLanczos
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_DRPLanczos.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_FGMRES
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_FGMRES.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_GMRESR
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_GMRESR.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_IPCG
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_IPCG.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_LBGMRESR
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_LBGMRESR.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_MINRES
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_MINRES.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_PCG
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_PCG.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_PLanczos
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_PLanczos.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

oops_add_test( TESTNAME simplifiedl95_RPCG
               MODELNAME l95
               YAMLNAME testinput/simplifiedl95_RPCG.yaml
               CTOL 1e-8
               EXENAME l95_4dvar.x )

#oops_add_test( TESTNAME simplifiedl95_RPLanczos
#               MODELNAME l95
#               YAMLNAME testinput/simplifiedl95_RPLanczos.yaml
#               CTOL 1e-8
#               EXENAME l95_4dvar.x )
#
