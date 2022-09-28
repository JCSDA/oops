list( APPEND qg_testinput
  testinput/3densvar.yaml
  testinput/3dvar.yaml
  testinput/3dvar_change_var.yaml
  testinput/3dvar_full_inverse.yaml
  testinput/3dvar_hybrid.yaml
  testinput/3dvar_hybrid_wo_jb_evaluation.yaml
  testinput/3dfgat.yaml
  testinput/4densvar.yaml
  testinput/4densvar_hybrid.yaml
  testinput/4dvar_dripcg.yaml
  testinput/4dvar_drpcg_lmp.yaml
  testinput/4dvar_drpfom.yaml
  testinput/4dvar_drplanczos.yaml
  testinput/4dvar_drplanczos_hybrid.yaml
  testinput/4dvar_forcing.yaml
  testinput/4dvar_htlm_qgtlm_dripcg.yaml
  testinput/4dvar_htlm_qgtlm_dripcg_pert_heat.yaml
  testinput/4dvar_ipcg.yaml
  testinput/4dvar_obs_biased.yaml
  testinput/4dvar_rpcg.yaml
  #testinput/4dvar_saddlepoint.yaml
  testinput/addincrement.yaml
  testinput/addincrement_scaled.yaml
  testinput/analytic_forecast.yaml
  testinput/convertincrement.yaml
  testinput/convertstate.yaml
  testinput/dfi.yaml
  testinput/diffstates.yaml
  testinput/dirac_cov.yaml
  testinput/dirac_hyb_field.yaml
  testinput/dirac_hyb_value.yaml
  testinput/dirac_loc_3d.yaml
  testinput/dirac_loc_4d.yaml
  testinput/dirac_no_loc.yaml
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
  testinput/eda_4dvar_1.yaml
  testinput/eda_4dvar_2.yaml
  testinput/eda_4dvar_3.yaml
  testinput/eda_4dvar_4.yaml
  testinput/eda_4dvar.yaml
  testinput/ens_forecast_1.yaml
  testinput/ens_forecast_2.yaml
  testinput/ens_forecast.yaml
  testinput/ens_hofx_1.yaml
  testinput/ens_hofx_2.yaml
  testinput/ens_hofx_3.yaml
  testinput/ens_hofx_4.yaml
  testinput/ens_hofx_5.yaml
  testinput/ens_hofx.yaml
  testinput/ens_recenter.yaml
  testinput/ens_variance.yaml
  testinput/ens_variance_inflation_field.yaml
  testinput/ens_variance_inflation_value.yaml
  testinput/error_covariance.yaml
  testinput/forecast.yaml
  testinput/forecast_control_htlm_pert_heat.yaml
  testinput/gen_ens_pert_B.yaml
  testinput/gen_ens_pert_B_HTLM.yaml
  testinput/gen_ens_pert_B_HTLM_pert_heat.yaml
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
  testinput/hybrid_linear_model_pert_heat.yaml
  testinput/increment.yaml
  testinput/letkf.yaml
  testinput/linear_model.yaml
  testinput/linear_obsoperator.yaml
  testinput/linear_variable_change.yaml
  testinput/localization.yaml
  testinput/locations.yaml
  testinput/make_obs_3d.yaml
  testinput/make_obs_4d_12h.yaml
  testinput/make_obs_4d_24h.yaml
  testinput/make_obs_4d_24h_pert_heat.yaml
  testinput/make_obs_4d_24h_filter_ordering.yaml
  testinput/make_obs_4d_biased.yaml
  testinput/model.yaml
  testinput/modelauxcovariance.yaml
  testinput/modelauxincrement.yaml
  testinput/obs_aux.yaml
  testinput/obs_iterator.yaml
  testinput/obs_localization.yaml
  testinput/obsoperator.yaml
  testinput/obsspace.yaml
  testinput/obsvector.yaml
  testinput/rtpp.yaml
  testinput/state.yaml
  testinput/static_b_init.yaml
  testinput/truth.yaml
  testinput/truth_pert_heat.yaml
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
  testoutput/4densvar.test
  testoutput/4densvar_hybrid.test
  testoutput/4dvar_dripcg.test
  testoutput/4dvar_htlm_qgtlm_dripcg.test
  testoutput/4dvar_htlm_qgtlm_dripcg_pert_heat.test
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
  testoutput/convertstate.test
  testoutput/addincrement.test
  testoutput/addincrement_scaled.test
  testoutput/dfi.test
  testoutput/diffstates.test
  testoutput/dirac_cov.test
  testoutput/dirac_hyb_field.test
  testoutput/dirac_hyb_value.test
  testoutput/dirac_loc_3d.test
  testoutput/dirac_loc_4d.test
  testoutput/dirac_no_loc.test
  testoutput/eda_3dfgat.test
  testoutput/eda_3dvar.test
  testoutput/eda_3dvar_block.test
  testoutput/eda_4dvar.test
  testoutput/ens_forecast.test
  testoutput/ens_hofx.test
  testoutput/ens_recenter.test
  testoutput/ens_variance.test
  testoutput/ens_variance_inflation_field.test
  testoutput/ens_variance_inflation_value.test
  testoutput/forecast.test
  testoutput/forecast_control_htlm_pert_heat.test
  testoutput/gen_ens_pert_B.test
  testoutput/gen_ens_pert_B_HTLM.test
  testoutput/gen_ens_pert_B_HTLM_pert_heat.test
  testoutput/hofx.test
  testoutput/hofx_tinterp.test
  testoutput/hofx3d.test
  testoutput/hybridgain_analysis.test
  testoutput/hybridgain_increment.test
  testoutput/letkf.test
  testoutput/make_obs_3d.test
  testoutput/make_obs_4d_12h.test
  testoutput/make_obs_4d_24h.test
  testoutput/make_obs_4d_24h_filter_ordering.test
  testoutput/make_obs_4d_24h_pert_heat.test
  testoutput/make_obs_4d_biased.test
  testoutput/rtpp.test
  testoutput/static_b_init.test
  testoutput/truth.test
  testoutput/truth_pert_heat.test
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

ecbuild_add_test( TARGET test_qg_truth_pert_heat
                  OMP 2
                  ARGS testinput/truth_pert_heat.yaml
                  COMMAND qg_forecast.x )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h
                  OMP 2
                  ARGS testinput/make_obs_4d_24h.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h_filter_ordering
                  OMP 2
                  ARGS testinput/make_obs_4d_24h_filter_ordering.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h_pert_heat
                  OMP 2
                  ARGS testinput/make_obs_4d_24h_pert_heat.yaml
                  COMMAND  qg_hofx.x
                  TEST_DEPENDS test_qg_truth_pert_heat )


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
                  TEST_DEPENDS test_qg_truth )

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
                  TEST_DEPENDS test_qg_truth test_qg_forecast_control_htlm test_qg_gen_ens_pert_B_HTLM )

ecbuild_add_test( TARGET  test_qg_hybrid_linear_model_pert_heat
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/hybrid_linear_model_pert_heat.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_truth_pert_heat test_qg_forecast_control_htlm_pert_heat test_qg_gen_ens_pert_B_HTLM_pert_heat )

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
                  ARGS "testinput/obsspace.yaml"
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
                  ARGS "testinput/obsspace.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/obs_aux.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux_increment
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/obsspace.yaml"
                  LIBS    qg
                  TEST_DEPENDS test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET  test_qg_obs_aux_covariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/obsspace.yaml"
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
ecbuild_add_test( TARGET test_qg_forecast_control_htlm_pert_heat
                  OMP 2
                  ARGS testinput/forecast_control_htlm_pert_heat.yaml
                  COMMAND  qg_forecast.x
                  TEST_DEPENDS test_qg_truth_pert_heat )


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

#####################################################################
# ensemble-related tests
#####################################################################

ecbuild_add_test( TARGET test_qg_gen_ens_pert_B
                  OMP 2
                  ARGS testinput/gen_ens_pert_B.yaml
                  COMMAND  qg_gen_ens_pert_B.x
                  TEST_DEPENDS test_qg_truth )
ecbuild_add_test( TARGET test_qg_gen_ens_pert_B_HTLM
                  OMP 2
                  ARGS testinput/gen_ens_pert_B_HTLM.yaml
                  COMMAND  qg_gen_ens_pert_B.x
                  TEST_DEPENDS test_qg_forecast test_qg_truth  )
                  
ecbuild_add_test( TARGET test_qg_gen_ens_pert_B_HTLM_pert_heat
                  OMP 2
                  ARGS testinput/gen_ens_pert_B_HTLM_pert_heat.yaml
                  COMMAND  qg_gen_ens_pert_B.x
                  TEST_DEPENDS test_qg_forecast_control_htlm_pert_heat test_qg_truth_pert_heat  )

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

ecbuild_add_test( TARGET test_qg_ens_variance
                  OMP 2
                  ARGS testinput/ens_variance.yaml
                  COMMAND  qg_ens_variance.x
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

ecbuild_add_test( TARGET test_qg_ens_variance_inflation_field
                  OMP 2
                  ARGS testinput/ens_variance_inflation_field.yaml
                  COMMAND  qg_ens_variance.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B test_qg_uniform_field_inflation )

ecbuild_add_test( TARGET test_qg_ens_variance_inflation_value
                  OMP 2
                  ARGS testinput/ens_variance_inflation_value.yaml
                  COMMAND  qg_ens_variance.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_rtpp
                  ARGS testinput/rtpp.yaml
                  COMMAND  qg_rtpp.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B )

#####################################################################
# other tests
#####################################################################

ecbuild_add_test( TARGET test_qg_static_b_init
                  MPI 1
                  OMP 2
                  ARGS testinput/static_b_init.yaml
                  COMMAND  qg_staticbinit.x
                  TEST_DEPENDS test_qg_truth )

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

#####################################################################
# QG dirac tests
#####################################################################

ecbuild_add_test( TARGET test_qg_dirac_cov
                  OMP 2
                  ARGS testinput/dirac_cov.yaml
                  COMMAND  qg_dirac.x
                  TEST_DEPENDS test_qg_forecast )

ecbuild_add_test( TARGET test_qg_dirac_hyb_value
                  OMP 2
                  ARGS testinput/dirac_hyb_value.yaml
                  COMMAND  qg_dirac.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_uniform_field_hybrid
                  OMP 2
                  ARGS testinput/uniform_field_hybrid.yaml
                  COMMAND qg_forecast.x )

ecbuild_add_test( TARGET test_qg_dirac_hyb_field
                  OMP 2
                  ARGS testinput/dirac_hyb_field.yaml
                  COMMAND  qg_dirac.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B test_qg_uniform_field_hybrid )

ecbuild_add_test( TARGET test_qg_dirac_loc_3d
                  OMP 2
                  ARGS testinput/dirac_loc_3d.yaml
                  COMMAND  qg_dirac.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_dirac_loc_4d
                  OMP 2
                  MPI 13
                  ARGS testinput/dirac_loc_4d.yaml
                  COMMAND  qg_dirac.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )

ecbuild_add_test( TARGET test_qg_dirac_no_loc
                  OMP 2
                  ARGS testinput/dirac_no_loc.yaml
                  COMMAND  qg_dirac.x
                  TEST_DEPENDS test_qg_forecast test_qg_gen_ens_pert_B )


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

ecbuild_add_test( TARGET test_qg_4dvar_htlm_qgtlm_dripcg
                  OMP 2
                  ARGS testinput/4dvar_htlm_qgtlm_dripcg.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_forecast test_qg_make_obs_4d_24h )

ecbuild_add_test( TARGET test_qg_4dvar_htlm_qgtlm_dripcg_pert_heat
                  OMP 2
                  ARGS testinput/4dvar_htlm_qgtlm_dripcg_pert_heat.yaml
                  COMMAND  qg_4dvar.x
                  TEST_DEPENDS test_qg_gen_ens_pert_B_HTLM_pert_heat test_qg_make_obs_4d_24h_pert_heat)

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

#####################################################################
# LETKF tests
#####################################################################

ecbuild_add_test( TARGET test_qg_letkf
                  ARGS testinput/letkf.yaml
                  COMMAND  qg_letkf.x
                  OMP 2
                  TEST_DEPENDS test_qg_make_obs_3d test_qg_gen_ens_pert_B )
