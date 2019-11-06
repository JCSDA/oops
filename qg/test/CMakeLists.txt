list( APPEND qg_test_input
  testinput/3densvar.yaml
  testinput/3densvar_bump.yaml
  testinput/3dvar.yaml
  testinput/3dvar_bump.yaml
  testinput/3dvar_change_var.yaml
  testinput/3dvar_hybrid.yaml
  testinput/3dvar_hybrid_bump.yaml
  testinput/3dfgat.yaml
  testinput/4densvar.yaml
  testinput/4densvar_bump.yaml
  testinput/4dvar_dripcg.yaml
  testinput/4dvar_drpcg_lmp.yaml
  testinput/4dvar_drpfom.yaml
  testinput/4dvar_drplanczos.yaml
  testinput/4dvar_drplanczos_bump.yaml
  testinput/4dvar_drplanczos_hybrid.yaml
  testinput/4dvar_drplanczos_hybrid_bump.yaml
  testinput/4dvar_forcing.yaml
  testinput/4dvar_ipcg.yaml
  testinput/4dvar_obs_biased.yaml
  testinput/4dvar_rpcg.yaml
  testinput/4dvar_saddlepoint.yaml
  testinput/analytic_forecast.yaml
  testinput/convertstate.yaml
  testinput/dfi.yaml
  testinput/dirac_bump_cov.yaml
  testinput/dirac_bump_hyb.yaml
  testinput/dirac_bump_lct.yaml
  testinput/dirac_bump_loc_3d.yaml
  testinput/dirac_bump_loc_4d.yaml
  testinput/dirac_qg_cov.yaml
  testinput/dirac_qg_hyb.yaml
  testinput/dirac_qg_loc_3d.yaml
  testinput/dirac_qg_loc_4d.yaml
  testinput/dirac_no_loc.yaml
  testinput/eda_4dvar_pert_obs.yaml
  testinput/ens_forecast.yaml
  testinput/ens_hofx.yaml
  testinput/ens_recenter.yaml
  testinput/ens_variance.yaml
  testinput/forecast.yaml
  testinput/gen_ens_pert_B.yaml
  testinput/hofx.yaml
  testinput/hofx3d.yaml
  testinput/interfaces.yaml
  testinput/letkf.yaml
  testinput/make_obs_3d.yaml
  testinput/make_obs_4d_12h.yaml
  testinput/make_obs_4d_24h.yaml
  testinput/make_obs_4d_biased.yaml
  testinput/parameters_bump_cov.yaml
  testinput/parameters_bump_hyb.yaml
  testinput/parameters_bump_lct.yaml
  testinput/parameters_bump_loc_3d.yaml
  testinput/parameters_bump_loc_4d.yaml
  testinput/static_b_init.yaml
  testinput/truth.yaml
)

list( APPEND qg_testoutput
  testoutput/3densvar.test
  testoutput/3densvar_bump.test
  testoutput/3dvar.test
  testoutput/3dvar_bump.test
  testoutput/3dvar_change_var.test
  testoutput/3dvar_hybrid.test
  testoutput/3dvar_hybrid_bump.test
  testoutput/3dfgat.test
  testoutput/4densvar.test
  testoutput/4densvar_bump.test
  testoutput/4dvar_dripcg.test
  testoutput/4dvar_drpcg_lmp.test
  testoutput/4dvar_drpfom.test
  testoutput/4dvar_drplanczos.test
  testoutput/4dvar_drplanczos_bump.test
  testoutput/4dvar_drplanczos_hybrid.test
  testoutput/4dvar_drplanczos_hybrid_bump.test
  testoutput/4dvar_forcing.test
  testoutput/4dvar_ipcg.test
  testoutput/4dvar_obs_biased.test
  testoutput/4dvar_rpcg.test
  testoutput/4dvar_saddlepoint.test
  testoutput/analytic_forecast.test
  testoutput/dfi.test
  testoutput/dirac_bump_cov.test
  testoutput/dirac_bump_hyb.test
  testoutput/dirac_bump_lct.test
  testoutput/dirac_bump_loc_3d.test
  testoutput/dirac_bump_loc_4d.test
  testoutput/dirac_qg_cov.test
  testoutput/dirac_qg_hyb.test
  testoutput/dirac_qg_loc_3d.test
  testoutput/dirac_qg_loc_4d.test
  testoutput/dirac_no_loc.test
#  testoutput/eda_4dvar_pert_obs.test
  testoutput/ens_forecast.test
  testoutput/ens_hofx.test
  testoutput/ens_recenter.test
  testoutput/ens_variance.test
  testoutput/forecast.test
  testoutput/gen_ens_pert_B.test
  testoutput/hofx.test
  testoutput/hofx3d.test
  testoutput/letkf.test
  testoutput/make_obs_3d.test
  testoutput/make_obs_4d_12h.test
  testoutput/make_obs_4d_24h.test
  testoutput/make_obs_4d_biased.test
  testoutput/parameters_bump_cov.test
  testoutput/parameters_bump_hyb.test
  testoutput/parameters_bump_lct.test
  testoutput/parameters_bump_loc_3d.test
  testoutput/parameters_bump_loc_4d.test
  testoutput/static_b_init.test
  testoutput/truth.test
)

# Create data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${qg_test_input})
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
                       ${qg_test_input}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)

#####################################################################
# truth and make_obs_4d_24h are required for interface tests
#####################################################################

ecbuild_add_test( TARGET test_qg_truth_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_forecast.x
                  ARGS testinput/truth.yaml testoutput/truth.log.out
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_truth_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/truth.log.out testoutput/truth.test
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_make_obs.x
                  ARGS testinput/make_obs_4d_24h.yaml testoutput/make_obs_4d_24h.log.out
                  DEPENDS qg_make_obs.x
                  TEST_DEPENDS test_qg_truth_run )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/make_obs_4d_24h.log.out testoutput/make_obs_4d_24h.test
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#####################################################################
# interface tests
#####################################################################

ecbuild_add_test( TARGET  test_qg_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_geometry_iterator
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_state
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_model_aux
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_model
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_error_covariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_modelauxincrement
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run
                  ENABLED OFF)

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_model_aux_covariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_linear_model
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_locations
                  SOURCES executables/TestLocations.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_obsspace
                  SOURCES executables/TestObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_localobsspace
                  SOURCES executables/TestLocalObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_obs_vector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_obs_operator
                  SOURCES executables/TestObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_lin_obs_operator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_obs_error_covar
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_obs_aux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_obs_aux_increment
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_obs_aux_covariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET  test_qg_linear_variable_change
                  MPI    1
                  SOURCES executables/TestLinearVariableChange.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_run )

#####################################################################
# forecast-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_analytic_forecast_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_forecast.x
                  ARGS testinput/analytic_forecast.yaml testoutput/analytic_forecast.log.out
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_analytic_forecast_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/analytic_forecast.log.out testoutput/analytic_forecast.test
                  TEST_DEPENDS test_qg_analytic_forecast_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_forecast_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_forecast.x
                  ARGS testinput/forecast.yaml testoutput/forecast.log.out
                  DEPENDS qg_forecast.x
                  TEST_DEPENDS test_qg_truth_run )

ecbuild_add_test( TARGET test_qg_forecast_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/forecast.log.out testoutput/forecast.test
                  TEST_DEPENDS test_qg_forecast_run )

#####################################################################
# obs-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_make_obs_3d_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_make_obs.x
                  ARGS testinput/make_obs_3d.yaml testoutput/make_obs_3d.log.out
                  DEPENDS qg_make_obs.x
                  TEST_DEPENDS test_qg_truth_run )

ecbuild_add_test( TARGET test_qg_make_obs_3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/make_obs_3d.log.out testoutput/make_obs_3d.test
                  TEST_DEPENDS test_qg_make_obs_3d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_make_obs_4d_12h_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_make_obs.x
                  ARGS testinput/make_obs_4d_12h.yaml testoutput/make_obs_4d_12h.log.out
                  DEPENDS qg_make_obs.x
                  TEST_DEPENDS test_qg_truth_run )

ecbuild_add_test( TARGET test_qg_make_obs_4d_12h_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/make_obs_4d_12h.log.out testoutput/make_obs_4d_12h.test
                  TEST_DEPENDS test_qg_make_obs_4d_12h_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_make_obs_4d_biased_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_make_obs.x
                  ARGS testinput/make_obs_4d_biased.yaml testoutput/make_obs_4d_biased.log.out
                  DEPENDS qg_make_obs.x
                  TEST_DEPENDS test_qg_truth_run )

ecbuild_add_test( TARGET test_qg_make_obs_4d_biased_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/make_obs_4d_biased.log.out testoutput/make_obs_4d_biased.test
                  TEST_DEPENDS test_qg_make_obs_4d_biased_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_hofx_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_hofx.x
                  ARGS testinput/hofx.yaml testoutput/hofx.log.out
                  DEPENDS qg_hofx.x
                  TEST_DEPENDS test_qg_make_obs_4d_12h_run )

ecbuild_add_test( TARGET test_qg_hofx_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/hofx.log.out testoutput/hofx.test
                  TEST_DEPENDS test_qg_hofx_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_hofx3d_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_hofx3d.x
                  ARGS testinput/hofx3d.yaml testoutput/hofx3d.log.out
                  DEPENDS qg_hofx3d.x
                  TEST_DEPENDS test_qg_make_obs_4d_12h_run )

ecbuild_add_test( TARGET test_qg_hofx3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/hofx3d.log.out testoutput/hofx3d.test
                  TEST_DEPENDS test_qg_hofx3d_run )

#####################################################################
# ensemble-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_gen_ens_pert_B_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_gen_ens_pert_B.x
                  ARGS testinput/gen_ens_pert_B.yaml testoutput/gen_ens_pert_B.log.out
                  DEPENDS qg_gen_ens_pert_B.x
                  TEST_DEPENDS test_qg_truth_run )

ecbuild_add_test( TARGET test_qg_gen_ens_pert_B_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/gen_ens_pert_B.log.out testoutput/gen_ens_pert_B.test
                  TEST_DEPENDS test_qg_gen_ens_pert_B_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_ens_forecast_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_ens_forecast.x
                  ARGS testinput/ens_forecast.yaml testoutput/ens_forecast.log.out
                  DEPENDS qg_ens_forecast.x
                  TEST_DEPENDS test_qg_forecast_run )

ecbuild_add_test( TARGET test_qg_ens_forecast_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/ens_forecast.log.out testoutput/ens_forecast.test
                  TEST_DEPENDS test_qg_ens_forecast_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_ens_hofx_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_ens_hofx.x
                  ARGS testinput/ens_hofx.yaml testoutput/ens_hofx.log.out
                  DEPENDS qg_ens_hofx.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B_run test_qg_make_obs_4d_12h_run )

ecbuild_add_test( TARGET test_qg_ens_hofx_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/ens_hofx.log.out testoutput/ens_hofx.test
                  TEST_DEPENDS test_qg_ens_hofx_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_ens_variance_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_ens_variance.x
                  ARGS testinput/ens_variance.yaml testoutput/ens_variance.log.out
                  DEPENDS qg_ens_variance.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_ens_variance_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/ens_variance.log.out testoutput/ens_variance.test
                  TEST_DEPENDS test_qg_ens_variance_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_ens_recenter_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_ens_recenter.x
                  ARGS testinput/ens_recenter.yaml testoutput/ens_recenter.log.out
                  DEPENDS qg_ens_recenter.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_ens_recenter_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/ens_recenter.log.out testoutput/ens_recenter.test
                  TEST_DEPENDS test_qg_ens_recenter_run )

#####################################################################
# other tests
#####################################################################

ecbuild_add_test( TARGET test_qg_static_b_init_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_staticbinit.x
                  ARGS testinput/static_b_init.yaml testoutput/static_b_init.log.out
                  DEPENDS qg_staticbinit.x
                  TEST_DEPENDS test_qg_truth_run )

ecbuild_add_test( TARGET test_qg_static_b_init_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/static_b_init.log.out testoutput/static_b_init.test
                  TEST_DEPENDS test_qg_static_b_init_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dfi_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dfi.x
                  ARGS testinput/dfi.yaml testoutput/dfi.log.out
                  DEPENDS qg_dfi.x
                  TEST_DEPENDS test_qg_forecast_run )

ecbuild_add_test( TARGET test_qg_dfi_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dfi.log.out testoutput/dfi.test
                  TEST_DEPENDS test_qg_dfi_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_convertstate_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_convertstate.x
                  ARGS testinput/convertstate.yaml testoutput/convertstate.log.out
                  DEPENDS qg_convertstate.x
                  TEST_DEPENDS test_qg_forecast_run )

#####################################################################
# BUMP parameters tests
#####################################################################

ecbuild_add_test( TARGET test_qg_parameters_bump_cov_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x
                  ARGS testinput/parameters_bump_cov.yaml testoutput/parameters_bump_cov.log.out
                  DEPENDS qg_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_parameters_bump_cov_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_cov.log.out testoutput/parameters_bump_cov.test
                  TEST_DEPENDS test_qg_parameters_bump_cov_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_parameters_bump_hyb_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x
                  ARGS testinput/parameters_bump_hyb.yaml testoutput/parameters_bump_hyb.log.out
                  DEPENDS qg_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_parameters_bump_hyb_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_hyb.log.out testoutput/parameters_bump_hyb.test
                  TEST_DEPENDS test_qg_parameters_bump_hyb_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_parameters_bump_lct_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x
                  ARGS testinput/parameters_bump_lct.yaml testoutput/parameters_bump_lct.log.out
                  DEPENDS qg_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_parameters_bump_lct_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_lct.log.out testoutput/parameters_bump_lct.test
                  TEST_DEPENDS test_qg_parameters_bump_lct_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_parameters_bump_loc_3d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x
                  ARGS testinput/parameters_bump_loc_3d.yaml testoutput/parameters_bump_loc_3d.log.out
                  DEPENDS qg_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_parameters_bump_loc_3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_loc_3d.log.out testoutput/parameters_bump_loc_3d.test
                  TEST_DEPENDS test_qg_parameters_bump_loc_3d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_parameters_bump_loc_4d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x
                  ARGS testinput/parameters_bump_loc_4d.yaml testoutput/parameters_bump_loc_4d.log.out
                  DEPENDS qg_estimate_parameters.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_parameters_bump_loc_4d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/parameters_bump_loc_4d.log.out testoutput/parameters_bump_loc_4d.test
                  TEST_DEPENDS test_qg_parameters_bump_loc_4d_run )

#####################################################################
# BUMP links
#####################################################################

ecbuild_add_test( TARGET test_qg_bump_links
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/bump_links.sh
                  TEST_DEPENDS test_qg_parameters_bump_cov_run
                               test_qg_parameters_bump_hyb_run
                               test_qg_parameters_bump_lct_run
                               test_qg_parameters_bump_loc_3d_run
                               test_qg_parameters_bump_loc_4d_run )

#####################################################################
# BUMP dirac tests
#####################################################################

ecbuild_add_test( TARGET test_qg_dirac_bump_cov_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_cov.yaml testoutput/dirac_bump_cov.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_dirac_bump_cov_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_cov.log.out testoutput/dirac_bump_cov.test
                  TEST_DEPENDS test_qg_dirac_bump_cov_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_bump_hyb_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_hyb.yaml testoutput/dirac_bump_hyb.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_dirac_bump_hyb_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_hyb.log.out testoutput/dirac_bump_hyb.test
                  TEST_DEPENDS test_qg_dirac_bump_hyb_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_bump_lct_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_lct.yaml testoutput/dirac_bump_lct.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_dirac_bump_lct_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_lct.log.out testoutput/dirac_bump_lct.test
                  TEST_DEPENDS test_qg_dirac_bump_lct_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_bump_loc_3d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_loc_3d.yaml testoutput/dirac_bump_loc_3d.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_dirac_bump_loc_3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_loc_3d.log.out testoutput/dirac_bump_loc_3d.test
                  TEST_DEPENDS test_qg_dirac_bump_loc_3d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_bump_loc_4d_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_bump_loc_4d.yaml testoutput/dirac_bump_loc_4d.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_dirac_bump_loc_4d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_bump_loc_4d.log.out testoutput/dirac_bump_loc_4d.test
                  TEST_DEPENDS test_qg_dirac_bump_loc_4d_run )

#####################################################################
# QG dirac tests
#####################################################################

ecbuild_add_test( TARGET test_qg_dirac_qg_cov_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_qg_cov.yaml testoutput/dirac_qg_cov.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run )

ecbuild_add_test( TARGET test_qg_dirac_qg_cov_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_qg_cov.log.out testoutput/dirac_qg_cov.test
                  TEST_DEPENDS test_qg_dirac_qg_cov_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_qg_hyb_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_qg_hyb.yaml testoutput/dirac_qg_hyb.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_dirac_qg_hyb_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_qg_hyb.log.out testoutput/dirac_qg_hyb.test
                  TEST_DEPENDS test_qg_dirac_qg_hyb_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_qg_loc_3d_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_qg_loc_3d.yaml testoutput/dirac_qg_loc_3d.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_dirac_qg_loc_3d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_qg_loc_3d.log.out testoutput/dirac_qg_loc_3d.test
                  TEST_DEPENDS test_qg_dirac_qg_loc_3d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_qg_loc_4d_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_qg_loc_4d.yaml testoutput/dirac_qg_loc_4d.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_dirac_qg_loc_4d_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_qg_loc_4d.log.out testoutput/dirac_qg_loc_4d.test
                  TEST_DEPENDS test_qg_dirac_qg_loc_4d_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_dirac_no_loc_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_dirac.x
                  ARGS testinput/dirac_no_loc.yaml testoutput/dirac_no_loc.log.out
                  DEPENDS qg_dirac.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run )

ecbuild_add_test( TARGET test_qg_dirac_no_loc_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/dirac_no_loc.log.out testoutput/dirac_no_loc.test
                  TEST_DEPENDS test_qg_dirac_no_loc_run )

#####################################################################
# 3d variational tests
#####################################################################

ecbuild_add_test( TARGET test_qg_3densvar_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3densvar.yaml testoutput/3densvar.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_3d_run )

ecbuild_add_test( TARGET test_qg_3densvar_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3densvar.log.out testoutput/3densvar.test
                  TEST_DEPENDS test_qg_3densvar_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3densvar_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3densvar_bump.yaml testoutput/3densvar_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_3d_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_3densvar_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3densvar_bump.log.out testoutput/3densvar_bump.test
                  TEST_DEPENDS test_qg_3densvar_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3dvar_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dvar.yaml testoutput/3dvar.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_3d_run )

ecbuild_add_test( TARGET test_qg_3dvar_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar.log.out testoutput/3dvar.test
                  TEST_DEPENDS test_qg_3dvar_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3dvar_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dvar_bump.yaml testoutput/3dvar_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_3d_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_3dvar_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar_bump.log.out testoutput/3dvar_bump.test
                  TEST_DEPENDS test_qg_3dvar_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3dvar_change_var_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dvar_change_var.yaml testoutput/3dvar_change_var.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_3d_run )

ecbuild_add_test( TARGET test_qg_3dvar_change_var_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar_change_var.log.out testoutput/3dvar_change_var.test
                  TEST_DEPENDS test_qg_3dvar_change_var_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3dvar_hybrid_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dvar_hybrid.yaml testoutput/3dvar_hybrid.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_3d_run )

ecbuild_add_test( TARGET test_qg_3dvar_hybrid_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar_hybrid.log.out testoutput/3dvar_hybrid.test
                  TEST_DEPENDS test_qg_3dvar_hybrid_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3dvar_hybrid_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dvar_hybrid_bump.yaml testoutput/3dvar_hybrid_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_3d_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_3dvar_hybrid_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dvar_hybrid_bump.log.out testoutput/3dvar_hybrid_bump.test
                  TEST_DEPENDS test_qg_3dvar_hybrid_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3dfgat_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/3dfgat.yaml testoutput/3dfgat.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_3d_run )

ecbuild_add_test( TARGET test_qg_3dfgat_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/3dfgat.log.out testoutput/3dfgat.test
                  TEST_DEPENDS test_qg_3dfgat_run )

#####################################################################
# 4d variational tests
#####################################################################

ecbuild_add_test( TARGET test_qg_4densvar_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4densvar.yaml testoutput/4densvar.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_4d_12h_run )

ecbuild_add_test( TARGET test_qg_4densvar_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4densvar.log.out testoutput/4densvar.test
                  TEST_DEPENDS test_qg_4densvar_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4densvar_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4densvar_bump.yaml testoutput/4densvar_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_4d_12h_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_4densvar_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4densvar_bump.log.out testoutput/4densvar_bump.test
                  TEST_DEPENDS test_qg_4densvar_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_dripcg_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_dripcg.yaml testoutput/4dvar_dripcg.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_dripcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_dripcg.log.out testoutput/4dvar_dripcg.test
                  TEST_DEPENDS test_qg_4dvar_dripcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_drpcg_lmp_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drpcg_lmp.yaml testoutput/4dvar_drpcg_lmp.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_drpcg_lmp_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drpcg_lmp.log.out testoutput/4dvar_drpcg_lmp.test
                  TEST_DEPENDS test_qg_4dvar_drpcg_lmp_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_drpfom_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drpfom.yaml testoutput/4dvar_drpfom.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_drpfom_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drpfom.log.out testoutput/4dvar_drpfom.test
                  TEST_DEPENDS test_qg_4dvar_drpfom_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drplanczos.yaml testoutput/4dvar_drplanczos.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drplanczos.log.out testoutput/4dvar_drplanczos.test
                  TEST_DEPENDS test_qg_4dvar_drplanczos_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drplanczos_bump.yaml testoutput/4dvar_drplanczos_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_4d_24h_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drplanczos_bump.log.out testoutput/4dvar_drplanczos_bump.test
                  TEST_DEPENDS test_qg_4dvar_drplanczos_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drplanczos_hybrid.yaml testoutput/4dvar_drplanczos_hybrid.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drplanczos_hybrid.log.out testoutput/4dvar_drplanczos_hybrid.test
                  TEST_DEPENDS test_qg_4dvar_drplanczos_hybrid_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid_bump_run
                  MPI    1
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_drplanczos_hybrid_bump.yaml testoutput/4dvar_drplanczos_hybrid_bump.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_gen_ens_pert_B_run test_qg_make_obs_4d_24h_run test_qg_bump_links )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid_bump_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_drplanczos_hybrid_bump.log.out testoutput/4dvar_drplanczos_hybrid_bump.test
                  TEST_DEPENDS test_qg_4dvar_drplanczos_hybrid_bump_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_forcing_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_forcing.yaml testoutput/4dvar_forcing.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_forcing_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_forcing.log.out testoutput/4dvar_forcing.test
                  TEST_DEPENDS test_qg_4dvar_forcing_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_ipcg_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_ipcg.yaml testoutput/4dvar_ipcg.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_ipcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_ipcg.log.out testoutput/4dvar_ipcg.test
                  TEST_DEPENDS test_qg_4dvar_ipcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_obs_biased_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_obs_biased.yaml testoutput/4dvar_obs_biased.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run  test_qg_make_obs_4d_biased_run )

ecbuild_add_test( TARGET test_qg_4dvar_obs_biased_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_obs_biased.log.out testoutput/4dvar_obs_biased.test
                  TEST_DEPENDS test_qg_4dvar_obs_biased_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_rpcg_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_rpcg.yaml testoutput/4dvar_rpcg.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run  test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_rpcg_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_rpcg.log.out testoutput/4dvar_rpcg.test
                  TEST_DEPENDS test_qg_4dvar_rpcg_run )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_4dvar_saddlepoint_run
                  OMP    2
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_4dvar.x
                  ARGS testinput/4dvar_saddlepoint.yaml testoutput/4dvar_saddlepoint.log.out
                  DEPENDS qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast_run  test_qg_make_obs_4d_24h_run )

ecbuild_add_test( TARGET test_qg_4dvar_saddlepoint_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/4dvar_saddlepoint.log.out testoutput/4dvar_saddlepoint.test
                  TEST_DEPENDS test_qg_4dvar_saddlepoint_run )

#####################################################################
# EDA tests
#####################################################################

ecbuild_add_test( TARGET test_qg_eda_4dvar_pert_obs_run
                  MPI    4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_4dvar_pert_obs.yaml
                  DEPENDS qg_eda.x
                  TEST_DEPENDS test_qg_forecast_run  test_qg_make_obs_4d_24h_run )

#####################################################################
# LETKF tests
#####################################################################

ecbuild_add_test( TARGET test_qg_letkf_run
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_letkf.x
                  ARGS testinput/letkf.yaml testoutput/letkf.log.out
          DEPENDS qg_letkf.x
          TEST_DEPENDS test_qg_make_obs_3d_run test_qg_gen_ens_pert_B_run)

ecbuild_add_test( TARGET test_qg_letkf_compare
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS testoutput/letkf.log.out testoutput/letkf.test
                  TEST_DEPENDS test_qg_letkf_run )
