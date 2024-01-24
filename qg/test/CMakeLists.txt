list( APPEND qg_testinput
  testinput/3densvar.yaml
  testinput/3dvar.yaml
  testinput/3dvar_change_var.yaml
  testinput/3dvar_full_inverse.yaml
  testinput/3dvar_hybrid.yaml
  testinput/3dvar_hybrid_wo_jb_evaluation.yaml
  testinput/3dfgat.yaml
  testinput/4dfgat.yaml
  testinput/4densvar.yaml
  testinput/4densvar_hybrid.yaml
  testinput/4dvar_dripcg.yaml
  testinput/4dvar_drpcg_lmp.yaml
  testinput/4dvar_drpfom.yaml
  testinput/4dvar_drplanczos.yaml
  testinput/4dvar_drplanczos_hybrid.yaml
  testinput/4dvar_forcing.yaml
  testinput/4dvar_hybrid_linear_model.yaml
  testinput/4dvar_ipcg.yaml
  testinput/4dvar_obs_biased.yaml
  testinput/4dvar_rpcg.yaml
  #testinput/4dvar_saddlepoint.yaml
  testinput/addincrement.yaml
  testinput/addincrement_scaled.yaml
  testinput/analytic_forecast.yaml
  testinput/convertincrement.yaml
  testinput/convertincrement_inv.yaml
  testinput/convertstate.yaml
  testinput/dfi.yaml
  testinput/diffstates.yaml
  testinput/3dvar_single-obs_cov.yaml
  testinput/3dvar_single-obs_hyb_field.yaml
  testinput/3dvar_single-obs_hyb_value.yaml
  testinput/3dvar_single-obs_loc_3d.yaml
  testinput/3dvar_single-obs_loc_3d_ensemble_transform.yaml
  testinput/3dvar_single-obs_loc_3d_inflation_field.yaml
  testinput/3dvar_single-obs_loc_3d_inflation_value.yaml
  testinput/4densvar_single-obs_loc_4d_fast.yaml
  testinput/4densvar_single-obs_loc_4d_aggressive.yaml
  testinput/4densvar_single-obs_loc_4d_standard.yaml
  testinput/4densvar_single-obs_loc_4d_time_decay_fast.yaml
  testinput/4densvar_single-obs_loc_4d_time_decay_standard.yaml
  testinput/4densvar_single-obs_no_loc.yaml
  testinput/eda_3dfgat_1.yaml
  testinput/eda_3dfgat_2.yaml
  testinput/eda_3dfgat_3.yaml
  testinput/eda_3dfgat_4.yaml
  testinput/eda_3dfgat.yaml
  testinput/eda_3dvar_2.yaml
  testinput/eda_3dvar.yaml
  testinput/eda_3dvar_block_1.yaml
  testinput/eda_3dvar_block_2.yaml
  testinput/eda_3dvar_block_3.yaml
  testinput/eda_3dvar_block_4.yaml
  testinput/eda_3dvar_block.yaml
  testinput/eda_3dvar_control_pert.yaml
  testinput/eda_4dvar_1.yaml
  testinput/eda_4dvar_2.yaml
  testinput/eda_4dvar_3.yaml
  testinput/eda_4dvar_4.yaml
  testinput/eda_4dvar.yaml
  testinput/eda_4dvar_control_pert.yaml
  testinput/ens_forecast_1.yaml
  testinput/ens_forecast_2.yaml
  testinput/ens_forecast.yaml
  testinput/ens_hofx_1.yaml
  testinput/ens_hofx_2.yaml
  testinput/ens_hofx_3.yaml
  testinput/ens_hofx_4.yaml
  testinput/ens_hofx_5.yaml
  testinput/ens_hofx.yaml
  testinput/ens_infl_inc.yaml
  testinput/ens_infl_states.yaml
  testinput/ens_mean_variance.yaml
  testinput/ens_recenter.yaml
  testinput/error_covariance.yaml
  testinput/fieldset4D_par.yaml
  testinput/fieldset4D_seq.yaml
  testinput/forecast.yaml
  testinput/forecast_5level_pert_heat.yaml
  testinput/gen_ens_pert_B.yaml
  testinput/gen_ens_pert_B_5level_pert_heat.yaml
  testinput/gen_hybrid_linear_model_coeffs_with_subtimestep.yaml
  testinput/geometry_iterator.yaml
  testinput/geometry.yaml
  testinput/geovals.yaml
  testinput/getvalues.yaml
  testinput/hofx.yaml
  testinput/hofx_tinterp.yaml
  testinput/hofx3d.yaml
  testinput/hybridgain_analysis.yaml
  testinput/hybridgain_increment.yaml
  testinput/hybrid_linear_model.yaml
  testinput/hybrid_linear_model_gen_ens.yaml
  testinput/increment.yaml
  testinput/letkf.yaml
  testinput/linear_model.yaml
  testinput/linear_obsoperator.yaml
  testinput/linear_variable_change.yaml
  testinput/linearization_error.yaml
  testinput/localization.yaml
  testinput/locations.yaml
  testinput/make_obs_3d.yaml
  testinput/make_obs_4d_12h.yaml
  testinput/make_obs_4d_24h.yaml
  testinput/make_obs_4d_24h_5level_pert_heat.yaml
  testinput/make_obs_4d_24h_filter_ordering.yaml
  testinput/make_obs_4d_biased.yaml
  testinput/make_single-obs_3d.yaml
  testinput/make_single-obs_4d_12h.yaml
  testinput/model.yaml
  testinput/modelauxcovariance.yaml
  testinput/modelauxincrement.yaml
  testinput/norm.yaml
  testinput/obs_aux.yaml
  testinput/obs_iterator.yaml
  testinput/obs_localization.yaml
  testinput/obsoperator.yaml
  testinput/obsdatavector.yaml
  testinput/obsspace.yaml
  testinput/obserror_covariance.yaml
  testinput/obs_aux_increment.yaml
  testinput/obs_aux_covariance.yaml
  testinput/obsvector.yaml
  testinput/state.yaml
  testinput/truth.yaml
  testinput/truth_5level_pert_heat.yaml
  testinput/variable_change.yaml
  testinput/verticallocev.yaml
  testinput/verticallocev_io.yaml
  testinput/uniform_field_hybrid.yaml
  testinput/uniform_field_inflation.yaml
)

list( APPEND qg_testoutput
  testoutput/3densvar.test
  testoutput/3dvar.test
  testoutput/3dvar_change_var.test
  testoutput/3dvar_full_inverse.test
  testoutput/3dvar_hybrid.test
  testoutput/3dvar_hybrid_wo_jb_evaluation.test
  testoutput/3dfgat.test
  testoutput/4dfgat.test
  testoutput/4densvar.test
  testoutput/4densvar_hybrid.test
  testoutput/4dvar_dripcg.test
  testoutput/4dvar_hybrid_linear_model.test
  testoutput/4dvar_drpcg_lmp.test
  testoutput/4dvar_drpfom.test
  testoutput/4dvar_drplanczos.test
  testoutput/4dvar_drplanczos_hybrid.test
  testoutput/4dvar_forcing.test
  testoutput/4dvar_ipcg.test
  testoutput/4dvar_obs_biased.test
  testoutput/4dvar_rpcg.test
  #testoutput/4dvar_saddlepoint.test
  testoutput/analytic_forecast.test
  testoutput/convertincrement.test
  testoutput/convertincrement_inv.test
  testoutput/convertstate.test
  testoutput/addincrement.test
  testoutput/addincrement_scaled.test
  testoutput/dfi.test
  testoutput/diffstates.test
  testoutput/3dvar_single-obs_cov.test
  testoutput/3dvar_single-obs_hyb_field.test
  testoutput/3dvar_single-obs_hyb_value.test
  testoutput/3dvar_single-obs_loc_3d.test
  testoutput/3dvar_single-obs_loc_3d_ensemble_transform.test
  testoutput/3dvar_single-obs_loc_3d_inflation_field.test
  testoutput/3dvar_single-obs_loc_3d_inflation_value.test
  testoutput/4densvar_single-obs_loc_4d.test
  testoutput/4densvar_single-obs_loc_4d_time_decay.test
  testoutput/4densvar_single-obs_no_loc.test
  testoutput/eda_3dfgat.test
  testoutput/eda_3dvar.test
  testoutput/eda_3dvar_block.test
  testoutput/eda_3dvar_control_pert.test
  testoutput/eda_4dvar.test
  testoutput/eda_4dvar_control_pert.test
  testoutput/ens_forecast.test
  testoutput/ens_hofx.test
  testoutput/ens_infl_inc.test
  testoutput/ens_infl_states.test
  testoutput/ens_mean_variance.test
  testoutput/ens_recenter.test
  testoutput/forecast.test
  testoutput/forecast_5level_pert_heat.test
  testoutput/gen_ens_pert_B.test
  testoutput/gen_ens_pert_B_5level_pert_heat.test
  testoutput/gen_hybrid_linear_model_coeffs_with_subtimestep.test
  testoutput/hofx.test
  testoutput/hofx_tinterp.test
  testoutput/hofx3d.test
  testoutput/hybridgain_analysis.test
  testoutput/hybridgain_increment.test
  testoutput/letkf.test
  testoutput/linearization_error.test
  testoutput/make_obs_3d.test
  testoutput/make_obs_4d_12h.test
  testoutput/make_obs_4d_24h.test
  testoutput/make_obs_4d_24h_5level_pert_heat.test
  testoutput/make_obs_4d_24h_filter_ordering.test
  testoutput/make_obs_4d_biased.test
  testoutput/make_single-obs_3d.test
  testoutput/make_single-obs_4d_12h.test
  testoutput/truth.test
  testoutput/truth_5level_pert_heat.test
  testoutput/uniform_field_hybrid.test
  testoutput/uniform_field_inflation.test
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

ecbuild_add_test( TARGET test_qg_truth
                  OMP 2
                  ARGS testinput/truth.yaml
                  COMMAND qg_forecast.x )

ecbuild_add_test( TARGET test_qg_truth_5level_pert_heat
                  OMP 2
                  ARGS testinput/truth_5level_pert_heat.yaml
                  COMMAND qg_forecast.x )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h
                  OMP 2
                  ARGS testinput/make_obs_4d_24h.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h_5level_pert_heat
                  OMP 2
                  ARGS testinput/make_obs_4d_24h_5level_pert_heat.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth_5level_pert_heat )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h_filter_ordering
                  OMP 2
                  ARGS testinput/make_obs_4d_24h_filter_ordering.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

#####################################################################
# interface tests
#####################################################################

ecbuild_add_test( TARGET  test_qg_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_geometry_iterator
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS    "testinput/geometry_iterator.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS    "testinput/geovals.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_state
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_model
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_error_covariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/error_covariance.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_verticallocev
                  SOURCES executables/TestVerticalLocEV.cc
                  ARGS    "testinput/verticallocev.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_verticallocev_io
                  SOURCES executables/TestVerticalLocEV.cc
                  ARGS    "testinput/verticallocev_io.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth test_qg_verticallocev )

ecbuild_add_test( TARGET  test_qg_modelauxincrement
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS    "testinput/modelauxincrement.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth
                  ENABLED OFF)

ecbuild_add_test( TARGET  test_qg_modelauxcovariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/modelauxcovariance.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_linear_model
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linear_model.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_hybrid_linear_model
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/hybrid_linear_model.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_gen_hybrid_linear_model_coeffs_with_subtimestep )

ecbuild_add_test( TARGET  test_qg_hybrid_linear_model_gen_ens
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/hybrid_linear_model_gen_ens.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_forecast_5level_pert_heat )

ecbuild_add_test( TARGET test_qg_locations
                  SOURCES executables/TestLocations.cc
                  ARGS "testinput/locations.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_obsspace
                  SOURCES executables/TestObsSpace.cc
                  ARGS "testinput/obsspace.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_obsvector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/obsvector.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_obsdatavector
                  SOURCES executables/TestObsDataVector.cc
                  ARGS "testinput/obsdatavector.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obsoperator
                  SOURCES executables/TestObsOperator.cc
                  ARGS    "testinput/obsoperator.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_linear_obsoperator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS    "testinput/linear_obsoperator.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_obserror_covariance
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/obserror_covariance.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/obs_aux.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux_increment
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/obs_aux_increment.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux_covariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/obs_aux_covariance.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_iterator
                  SOURCES executables/TestObsIterator.cc
                  ARGS    "testinput/obs_iterator.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_localization
                  SOURCES executables/TestObsLocalization.cc
                  ARGS    "testinput/obs_localization.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/localization.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_variable_change
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/variable_change.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_linear_variable_change
                  MPI    1
                  SOURCES executables/TestLinearVariableChange.cc
                  ARGS    "testinput/linear_variable_change.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET  test_qg_getvalues
                  SOURCES executables/TestGetValues.cc
                  ARGS    "testinput/getvalues.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth )

#####################################################################
# forecast-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_analytic_forecast
                  OMP 2
                  ARGS testinput/analytic_forecast.yaml
                  COMMAND  qg_forecast.x )

ecbuild_add_test( TARGET test_qg_forecast
                  OMP 2
                  ARGS testinput/forecast.yaml
                  COMMAND  qg_forecast.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_forecast_5level_pert_heat
                  OMP 2
                  ARGS testinput/forecast_5level_pert_heat.yaml
                  COMMAND  qg_forecast.x
                  TEST_DEPENDS test_qg_truth_5level_pert_heat )

#####################################################################
# obs-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_make_obs_3d
                  OMP 2
                  ARGS testinput/make_obs_3d.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_make_obs_4d_12h
                  OMP 2
                  ARGS testinput/make_obs_4d_12h.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_make_obs_4d_biased
                  OMP 2
                  ARGS testinput/make_obs_4d_biased.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_hofx
                  OMP 2
                  ARGS testinput/hofx.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_make_obs_4d_12h )

ecbuild_add_test( TARGET test_qg_hofx_tinterp
                  OMP 2
                  ARGS testinput/hofx_tinterp.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_make_obs_4d_12h )

ecbuild_add_test( TARGET test_qg_hofx3d
                  OMP 2
                  ARGS testinput/hofx3d.yaml
                  COMMAND  qg_hofx3d.x
                  TEST_DEPENDS test_qg_make_obs_4d_12h )

ecbuild_add_test( TARGET test_qg_make_single-obs_3d
                  OMP 2
                  ARGS testinput/make_single-obs_3d.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_make_single-obs_4d_12h
                  OMP 2
                  ARGS testinput/make_single-obs_4d_12h.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

#####################################################################
# 4D-related tests
#####################################################################

ecbuild_add_test( TARGET  test_qg_fieldset4D_seq
                  MPI     1
                  SOURCES executables/TestFieldSet4D.cc
                  ARGS    "testinput/fieldset4D_seq.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth test_qg_forecast )

ecbuild_add_test( TARGET  test_qg_fieldset4D_par
                  MPI     3
                  SOURCES executables/TestFieldSet4D.cc
                  ARGS    "testinput/fieldset4D_par.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth test_qg_forecast )

ecbuild_add_test( TARGET  test_qg_Norm
                  MPI     1
                  SOURCES executables/TestNorm.cc
                  ARGS    "testinput/norm.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth test_qg_forecast )

#####################################################################
# ensemble-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_gen_ens_pert_B
                  OMP 2
                  ARGS testinput/gen_ens_pert_B.yaml
                  COMMAND  qg_gen_ens_pert_B.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_gen_ens_pert_B_5level_pert_heat
                  OMP 2
                  ARGS testinput/gen_ens_pert_B_5level_pert_heat.yaml
                  COMMAND  qg_gen_ens_pert_B.x
                  TEST_DEPENDS test_qg_truth_5level_pert_heat )

ecbuild_add_test( TARGET test_qg_gen_hybrid_linear_model_coeffs_with_subtimestep
                  OMP 2
                  ARGS testinput/gen_hybrid_linear_model_coeffs_with_subtimestep.yaml
                  COMMAND  qg_gen_hybrid_linear_model_coeffs.x
                  TEST_DEPENDS test_qg_forecast_5level_pert_heat test_qg_gen_ens_pert_B_5level_pert_heat )

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

ecbuild_add_test( TARGET test_qg_ens_mean_variance
                  OMP 2
                  ARGS testinput/ens_mean_variance.yaml
                  COMMAND  qg_ens_mean_variance.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_ens_recenter
                  OMP 2
                  ARGS testinput/ens_recenter.yaml
                  COMMAND  qg_ens_recenter.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_hybridgain_analysis
                  OMP 2
                  ARGS testinput/hybridgain_analysis.yaml
                  COMMAND  qg_hybridgain.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_hybridgain_increment
                  OMP 2
                  ARGS testinput/hybridgain_increment.yaml
                  COMMAND  qg_hybridgain.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_uniform_field_inflation
                  OMP 2
                  ARGS testinput/uniform_field_inflation.yaml
                  COMMAND qg_forecast.x )

ecbuild_add_test( TARGET test_qg_ens_inflation_inc
                  ARGS testinput/ens_infl_inc.yaml
                  COMMAND  qg_ensemble_inflation.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B test_qg_eda_3dvar_control_pert)

ecbuild_add_test( TARGET test_qg_ens_inflation_states
                  ARGS testinput/ens_infl_states.yaml
                  COMMAND  qg_ensemble_inflation.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B test_qg_eda_4dvar)

#####################################################################
# other tests
#####################################################################

ecbuild_add_test( TARGET test_qg_dfi
                  OMP 2
                  ARGS testinput/dfi.yaml
                  COMMAND  qg_dfi.x
                  TEST_DEPENDS test_qg_forecast )

ecbuild_add_test( TARGET test_qg_convertstate
                  OMP    2
                  ARGS testinput/convertstate.yaml
                  COMMAND qg_convertstate.x
                  TEST_DEPENDS test_qg_forecast )

ecbuild_add_test( TARGET test_qg_linearization_error
                  OMP    2
                  ARGS testinput/linearization_error.yaml
                  COMMAND qg_linearization_error.x
                  TEST_DEPENDS test_qg_4dvar_drplanczos_hybrid )

#####################################################################
# QG single-obs tests
#####################################################################

ecbuild_add_test( TARGET test_qg_3dvar_single-obs_cov
                  OMP 2
                  ARGS testinput/3dvar_single-obs_cov.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_single-obs_3d )

ecbuild_add_test( TARGET test_qg_3dvar_single-obs_hyb_value
                  OMP 2
                  ARGS testinput/3dvar_single-obs_hyb_value.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_single-obs_3d test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_uniform_field_hybrid
                  OMP 2
                  ARGS testinput/uniform_field_hybrid.yaml
                  COMMAND qg_forecast.x )

ecbuild_add_test( TARGET test_qg_3dvar_single-obs_hyb_field
                  OMP 2
                  ARGS testinput/3dvar_single-obs_hyb_field.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_single-obs_3d test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_3dvar_single-obs_loc_3d
                  OMP 2
                  ARGS testinput/3dvar_single-obs_loc_3d.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_single-obs_3d test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_3dvar_single-obs_loc_3d_ensemble_transform
                  OMP 2
                  ARGS testinput/3dvar_single-obs_loc_3d_ensemble_transform.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_single-obs_3d test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_3dvar_single-obs_loc_3d_inflation_field
                  OMP 2
                  ARGS testinput/3dvar_single-obs_loc_3d_inflation_field.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_single-obs_3d test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_3dvar_single-obs_loc_3d_inflation_value
                  OMP 2
                  ARGS testinput/3dvar_single-obs_loc_3d_inflation_value.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_single-obs_3d test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_4densvar_single-obs_loc_4d_fast
                  MPI 13
                  ARGS testinput/4densvar_single-obs_loc_4d_fast.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_single-obs_4d_12h )

ecbuild_add_test( TARGET test_qg_4densvar_single-obs_loc_4d_aggressive
                  MPI 13
                  ARGS testinput/4densvar_single-obs_loc_4d_aggressive.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_single-obs_4d_12h )

ecbuild_add_test( TARGET test_qg_4densvar_single-obs_loc_4d_standard
                  MPI 13
                  ARGS testinput/4densvar_single-obs_loc_4d_standard.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_single-obs_4d_12h )

ecbuild_add_test( TARGET test_qg_4densvar_single-obs_loc_4d_time_decay_fast
                  MPI 13
                  ARGS testinput/4densvar_single-obs_loc_4d_time_decay_fast.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_single-obs_4d_12h )

ecbuild_add_test( TARGET test_qg_4densvar_single-obs_loc_4d_time_decay_standard
                  MPI 13
                  ARGS testinput/4densvar_single-obs_loc_4d_time_decay_standard.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_single-obs_4d_12h )

ecbuild_add_test( TARGET test_qg_4densvar_single-obs_no_loc
                  MPI 13
                  ARGS testinput/4densvar_single-obs_no_loc.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_single-obs_4d_12h )

#####################################################################
# 3d variational tests
#####################################################################

ecbuild_add_test( TARGET test_qg_3densvar
                  OMP 2
                  ARGS testinput/3densvar.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_3dvar
                  OMP 2
                  ARGS testinput/3dvar.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

#--------------------------------------------------------------------

ecbuild_add_test( TARGET test_qg_3dvar_change_var
                  OMP 2
                  ARGS testinput/3dvar_change_var.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_3dvar_full_inverse
                  OMP 2
                  ARGS testinput/3dvar_full_inverse.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_3dvar_hybrid
                  OMP 2
                  ARGS testinput/3dvar_hybrid.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_3dvar_hybrid_wo_jb_evaluation
                  OMP 2
                  ARGS testinput/3dvar_hybrid_wo_jb_evaluation.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_3dfgat
                  OMP 2
                  ARGS testinput/3dfgat.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_4dfgat
                  OMP 2
                  ARGS testinput/4dfgat.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

#####################################################################
# 4d variational tests
#####################################################################

ecbuild_add_test( TARGET test_qg_4densvar
                  MPI 7
                  ARGS testinput/4densvar.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_4d_12h )

ecbuild_add_test( TARGET test_qg_4densvar_hybrid
                  MPI 7
                  ARGS testinput/4densvar_hybrid.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_4d_12h )

ecbuild_add_test( TARGET test_qg_4dvar_dripcg
                  OMP 2
                  ARGS testinput/4dvar_dripcg.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_4dvar_drpcg_lmp
                  OMP 2
                  ARGS testinput/4dvar_drpcg_lmp.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_4dvar_drpfom
                  OMP 2
                  ARGS testinput/4dvar_drpfom.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos
                  OMP 2
                  ARGS testinput/4dvar_drplanczos.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid
                  OMP 2
                  ARGS testinput/4dvar_drplanczos_hybrid.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_make_obs_4d_24h )

#ecbuild_add_test( TARGET test_qg_4dvar_forcing
#                  OMP 2
#                  ARGS testinput/4dvar_forcing.yaml
#                  COMMAND  qg_4dvar.x
#                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_4dvar_hybrid_linear_model
                  OMP 2
                  ARGS testinput/4dvar_hybrid_linear_model.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_make_obs_4d_24h_5level_pert_heat test_qg_gen_hybrid_linear_model_coeffs_with_subtimestep )

ecbuild_add_test( TARGET test_qg_4dvar_ipcg
                  OMP 2
                  ARGS testinput/4dvar_ipcg.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_4dvar_obs_biased
                  OMP 2
                  ARGS testinput/4dvar_obs_biased.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_biased )

ecbuild_add_test( TARGET test_qg_4dvar_rpcg
                  OMP 2
                  ARGS testinput/4dvar_rpcg.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

#ecbuild_add_test( TARGET test_qg_4dvar_saddlepoint
#                  MPI 2
#                  ARGS testinput/4dvar_saddlepoint.yaml
#                  COMMAND  qg_4dvar.x
#                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

#####################################################################
# EDA tests
#####################################################################

ecbuild_add_test( TARGET test_qg_eda_3dfgat
                  MPI    4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_3dfgat.yaml
                  DEPENDS qg_eda.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_eda_3dvar
                  MPI    1
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_3dvar.yaml
                  DEPENDS qg_eda.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d )

ecbuild_add_test( TARGET test_qg_eda_4dvar
                  MPI    4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_4dvar.yaml
                  DEPENDS qg_eda.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_eda_3dvar_block
                  MPI    4
                  ARGS testinput/eda_3dvar_block.yaml
                  COMMAND qg_eda.x
                  TEST_DEPENDS test_qg_make_obs_3d test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_eda_3dvar_control_pert
                  MPI    5
                  ARGS testinput/eda_3dvar_control_pert.yaml
                  COMMAND qg_ControlPert.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_3d test_qg_gen_ens_pert_B)

ecbuild_add_test( TARGET test_qg_eda_4dvar_control_pert
                  MPI    4
                  ARGS testinput/eda_4dvar_control_pert.yaml
                  COMMAND qg_ControlPert.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h test_qg_gen_ens_pert_B)

#####################################################################
# state-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_diffstates
                  ARGS testinput/diffstates.yaml
                  COMMAND  qg_diffstates.x
                  TEST_DEPENDS test_qg_eda_3dvar test_qg_eda_4dvar
                  DEPENDS qg_eda.x )

ecbuild_add_test( TARGET test_qg_addincrement
                  ARGS testinput/addincrement.yaml
                  COMMAND  qg_addincrement.x
                  TEST_DEPENDS test_qg_diffstates )

ecbuild_add_test( TARGET test_qg_addincrement_scaled
                  ARGS testinput/addincrement_scaled.yaml
                  COMMAND  qg_addincrement.x
                  TEST_DEPENDS test_qg_diffstates )

ecbuild_add_test( TARGET test_qg_convertincrement
                  OMP    2
                  ARGS testinput/convertincrement.yaml
                  COMMAND qg_convertincrement.x
                  TEST_DEPENDS test_qg_diffstates )

ecbuild_add_test( TARGET test_qg_convertincrement_inv
                  OMP    2
                  ARGS testinput/convertincrement_inv.yaml
                  COMMAND qg_convertincrement.x
                  TEST_DEPENDS test_qg_diffstates )

#####################################################################
# LETKF tests
#####################################################################

ecbuild_add_test( TARGET test_qg_letkf
                  ARGS testinput/letkf.yaml
                  COMMAND  qg_letkf.x
                  OMP 2
                  TEST_DEPENDS test_qg_make_obs_3d test_qg_gen_ens_pert_B )
