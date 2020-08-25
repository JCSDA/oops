list( APPEND qg_testinput
  testinput/3densvar.yaml
  testinput/3dvar.yaml
  testinput/3dvar_change_var.yaml
  testinput/3dvar_hybrid.yaml
  testinput/3dfgat.yaml
  testinput/4densvar.yaml
  testinput/4dvar_dripcg.yaml
  testinput/4dvar_drpcg_lmp.yaml
  testinput/4dvar_drpfom.yaml
  testinput/4dvar_drplanczos.yaml
  testinput/4dvar_drplanczos_hybrid.yaml
  testinput/4dvar_forcing.yaml
  testinput/4dvar_ipcg.yaml
  testinput/4dvar_obs_biased.yaml
  testinput/4dvar_rpcg.yaml
  testinput/4dvar_saddlepoint.yaml
  testinput/analytic_forecast.yaml
  testinput/convertstate.yaml
  testinput/dfi.yaml
  testinput/diffstates.yaml
  testinput/dirac_cov.yaml
  testinput/dirac_hyb.yaml
  testinput/dirac_loc_3d.yaml
  testinput/dirac_loc_4d.yaml
  testinput/dirac_no_loc.yaml
  testinput/eda_3dfgat.yaml
  testinput/eda_3dvar_member2.yaml
  testinput/eda_4dvar.yaml
  testinput/ens_forecast.yaml
  testinput/ens_hofx.yaml
  testinput/ens_recenter.yaml
  testinput/ens_variance.yaml
  testinput/ens_variance_inflation_field.yaml
  testinput/ens_variance_inflation_value.yaml
  testinput/forecast.yaml
  testinput/gen_ens_pert_B.yaml
  testinput/getvalues.yaml
  testinput/hofx.yaml
  testinput/hofx3d.yaml
  testinput/interfaces.yaml
  testinput/letkf.yaml
  testinput/lineargetvalues.yaml
  testinput/make_obs_3d.yaml
  testinput/make_obs_4d_12h.yaml
  testinput/make_obs_4d_24h.yaml
  testinput/make_obs_4d_biased.yaml
  testinput/rtpp.yaml
  testinput/static_b_init.yaml
  testinput/truth.yaml
  testinput/uniform_field.yaml
)

list( APPEND qg_testoutput
  testoutput/3densvar.test
  testoutput/3dvar.test
  testoutput/3dvar_change_var.test
  testoutput/3dvar_hybrid.test
  testoutput/3dfgat.test
  testoutput/4densvar.test
  testoutput/4dvar_dripcg.test
  testoutput/4dvar_drpcg_lmp.test
  testoutput/4dvar_drpfom.test
  testoutput/4dvar_drplanczos.test
  testoutput/4dvar_drplanczos_hybrid.test
  testoutput/4dvar_forcing.test
  testoutput/4dvar_ipcg.test
  testoutput/4dvar_obs_biased.test
  testoutput/4dvar_rpcg.test
  testoutput/4dvar_saddlepoint.test
  testoutput/analytic_forecast.test
  testoutput/dfi.test
  testoutput/diffstates.test
  testoutput/dirac_cov.test
  testoutput/dirac_hyb.test
  testoutput/dirac_loc_3d.test
  testoutput/dirac_loc_4d.test
  testoutput/dirac_no_loc.test
#  testoutput/eda_4dvar.test
  testoutput/ens_forecast.test
  testoutput/ens_hofx.test
  testoutput/ens_recenter.test
  testoutput/ens_variance.test
  testoutput/ens_variance_inflation_field.test
  testoutput/ens_variance_inflation_value.test
  testoutput/forecast.test
  testoutput/gen_ens_pert_B.test
  testoutput/hofx.test
  testoutput/hofx3d.test
  testoutput/letkf.test
  testoutput/make_obs_3d.test
  testoutput/make_obs_4d_12h.test
  testoutput/make_obs_4d_24h.test
  testoutput/make_obs_4d_biased.test
  testoutput/rtpp.test
  testoutput/static_b_init.test
  testoutput/truth.test
  testoutput/uniform_field.test
)

# Create data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${qg_testinput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${qg_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   qg_test_scripts
                       SOURCES_PACK
                       ${qg_testinput}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)


#####################################################################
# truth and make_obs_4d_24h are required for interface tests
#####################################################################

oops_add_test( TESTNAME truth
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/truth.yaml
               EXENAME qg_forecast.x )

oops_add_test( TESTNAME make_obs_4d_24h
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/make_obs_4d_24h.yaml
               EXENAME  qg_hofx.x
               TEST_DEPENDS test_qg_truth )


#####################################################################
# interface tests
#####################################################################

ecbuild_add_test( TARGET  test_qg_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_geometry_iterator
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_state
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_model_aux
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_model
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_error_covariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_modelauxincrement
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth
                  ENABLED OFF)

ecbuild_add_test( TARGET  test_qg_model_aux_covariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_linear_model
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_locations
                  SOURCES executables/TestLocations.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_obsspace
                  SOURCES executables/TestObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_localobsspace
                  SOURCES executables/TestLocalObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_obs_vector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_departures_ensemble
                  SOURCES executables/TestDeparturesEnsemble.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_departures
                  SOURCES executables/TestDepartures.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_operator
                  SOURCES executables/TestObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_lin_obs_operator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_obs_error_covar
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux_increment
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux_covariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_linear_variable_change
                  MPI    1
                  SOURCES executables/TestLinearVariableChange.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_getvalues
                  SOURCES executables/TestGetValues.cc
                  ARGS    "testinput/getvalues.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_lineargetvalues
                  SOURCES executables/TestLinearGetValues.cc
                  ARGS    "testinput/lineargetvalues.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )


#####################################################################
# forecast-related tests
#####################################################################

oops_add_test( TESTNAME analytic_forecast
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/analytic_forecast.yaml
               EXENAME  qg_forecast.x )

oops_add_test( TESTNAME forecast
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/forecast.yaml
               EXENAME  qg_forecast.x
               TEST_DEPENDS test_qg_truth )


#####################################################################
# obs-related tests
#####################################################################

oops_add_test( TESTNAME make_obs_3d
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/make_obs_3d.yaml
               EXENAME  qg_hofx.x
               TEST_DEPENDS test_qg_truth )

oops_add_test( TESTNAME make_obs_4d_12h
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/make_obs_4d_12h.yaml
               EXENAME  qg_hofx.x
               TEST_DEPENDS test_qg_truth )

oops_add_test( TESTNAME make_obs_4d_biased
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/make_obs_4d_biased.yaml
               EXENAME  qg_hofx.x
               TEST_DEPENDS test_qg_truth )

oops_add_test( TESTNAME hofx
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/hofx.yaml
               EXENAME  qg_hofx.x
               TEST_DEPENDS test_qg_make_obs_4d_12h )

oops_add_test( TESTNAME hofx3d
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/hofx3d.yaml
               EXENAME  qg_hofx_nomodel.x
               TEST_DEPENDS test_qg_make_obs_4d_12h )


#####################################################################
# ensemble-related tests
#####################################################################

oops_add_test( TESTNAME gen_ens_pert_B
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/gen_ens_pert_B.yaml
               EXENAME  qg_gen_ens_pert_B.x
               TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_ens_forecast
                  MPI    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_ens_forecast.x
                  ARGS testinput/ens_forecast.yaml testoutput/ens_forecast.log.out
                  DEPENDS qg_ens_forecast.x
                  TEST_DEPENDS test_qg_forecast )

ecbuild_add_test( TARGET test_qg_ens_hofx
                  MPI    5
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_ens_hofx.x
                  ARGS testinput/ens_hofx.yaml testoutput/ens_hofx.log.out
                  DEPENDS qg_ens_hofx.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B test_qg_make_obs_4d_12h )

oops_add_test( TESTNAME ens_variance
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/ens_variance.yaml
               EXENAME  qg_ens_variance.x
               TEST_DEPENDS test_qg_gen_ens_pert_B )

oops_add_test( TESTNAME ens_recenter
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/ens_recenter.yaml
               EXENAME  qg_ens_recenter.x
               TEST_DEPENDS test_qg_gen_ens_pert_B )

oops_add_test( TESTNAME uniform_field
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/uniform_field.yaml
               EXENAME qg_forecast.x )

oops_add_test( TESTNAME ens_variance_inflation_field
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/ens_variance_inflation_field.yaml
               EXENAME  qg_ens_variance.x
               TEST_DEPENDS test_qg_gen_ens_pert_B test_qg_uniform_field )

oops_add_test( TESTNAME ens_variance_inflation_value
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/ens_variance_inflation_value.yaml
               EXENAME  qg_ens_variance.x
               TEST_DEPENDS test_qg_gen_ens_pert_B )

oops_add_test( TESTNAME rtpp
               MODELNAME qg
               YAMLNAME testinput/rtpp.yaml
               EXENAME  qg_rtpp.x
               TEST_DEPENDS test_qg_gen_ens_pert_B )

#####################################################################
# other tests
#####################################################################

oops_add_test( TESTNAME static_b_init
               MODELNAME qg
               MPI 1
               OMP 2
               YAMLNAME testinput/static_b_init.yaml
               EXENAME  qg_staticbinit.x
               TEST_DEPENDS test_qg_truth )

oops_add_test( TESTNAME dfi
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/dfi.yaml
               EXENAME  qg_dfi.x
               TEST_DEPENDS test_qg_forecast )

ecbuild_add_test( TARGET test_qg_convertstate
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_convertstate.x
                  ARGS testinput/convertstate.yaml testoutput/convertstate.log.out
                  DEPENDS qg_convertstate.x
                  TEST_DEPENDS test_qg_forecast )


#####################################################################
# QG dirac tests
#####################################################################

oops_add_test( TESTNAME dirac_cov
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/dirac_cov.yaml
               EXENAME  qg_dirac.x
               TEST_DEPENDS test_qg_forecast )

oops_add_test( TESTNAME dirac_hyb
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/dirac_hyb.yaml
               EXENAME  qg_dirac.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )

oops_add_test( TESTNAME dirac_loc_3d
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/dirac_loc_3d.yaml
               EXENAME  qg_dirac.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )

oops_add_test( TESTNAME dirac_loc_4d
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/dirac_loc_4d.yaml
               EXENAME  qg_dirac.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )

oops_add_test( TESTNAME dirac_no_loc
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/dirac_no_loc.yaml
               EXENAME  qg_dirac.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )


#####################################################################
# 3d variational tests
#####################################################################

oops_add_test( TESTNAME 3densvar
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/3densvar.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_3d )

oops_add_test( TESTNAME 3dvar
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/3dvar.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

# This test will be changed in next PR (May 2020)
# oops_add_test( TESTNAME 3dvar_change_var
#                MODELNAME qg
#                OMP 2
#                YAMLNAME testinput/3dvar_change_var.yaml
#                EXENAME  qg_4dvar.x
#                TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_3d
#                COMPARE)
# This test will be changed in next PR (May 2020)

#--------------------------------------------------------------------

oops_add_test( TESTNAME 3dvar_hybrid
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/3dvar_hybrid.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_3d )

oops_add_test( TESTNAME 3dfgat
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/3dfgat.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

#####################################################################
# 4d variational tests
#####################################################################

oops_add_test( TESTNAME 4densvar
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4densvar.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_4d_12h )

oops_add_test( TESTNAME 4dvar_dripcg
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_dripcg.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_drpcg_lmp
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_drpcg_lmp.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_drpfom
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_drpfom.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_drplanczos
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_drplanczos.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_drplanczos_hybrid
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_drplanczos_hybrid.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_forcing
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_forcing.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_ipcg
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_ipcg.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_obs_biased
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_obs_biased.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_biased )

oops_add_test( TESTNAME 4dvar_rpcg
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_rpcg.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

oops_add_test( TESTNAME 4dvar_saddlepoint
               MODELNAME qg
               OMP 2
               YAMLNAME testinput/4dvar_saddlepoint.yaml
               EXENAME  qg_4dvar.x
               TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h
               CTOL 0.0000000001 )


#####################################################################
# EDA tests
#####################################################################

ecbuild_add_test( TARGET test_qg_eda_3dfgat
                  MPI    4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_3dfgat.yaml
                  DEPENDS qg_eda.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_eda_3dvar_member2
                  MPI    1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_3dvar_member2.yaml
                  DEPENDS qg_eda.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_eda_4dvar
                  MPI    4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_4dvar.yaml
                  DEPENDS qg_eda.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )


#####################################################################
# state-related tests
#####################################################################

oops_add_test( TESTNAME diffstates
               MODELNAME qg
               YAMLNAME testinput/diffstates.yaml
               EXENAME  qg_diffstates.x
               TEST_DEPENDS test_qg_eda_3dvar_member2  test_qg_eda_4dvar 
               DEPENDS qg_eda.x )


#####################################################################
# LETKF tests
#####################################################################

oops_add_test( TESTNAME letkf
               MODELNAME qg
               YAMLNAME testinput/letkf.yaml
               EXENAME  qg_letkf.x
               TEST_DEPENDS test_qg_make_obs_3d test_qg_gen_ens_pert_B )
