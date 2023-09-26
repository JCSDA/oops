list( APPEND l95_test_input
  testinput/3dvar.yaml
  testinput/3dvar_4fsoi_dripcg.yaml
  testinput/3dvar_4fsoi_pcg.yaml
  testinput/3dvar_noobs.yaml
  testinput/3dvar_qc.yaml
  testinput/3dvar_qc_filter_ordering_pre.yaml
  testinput/3dvar_qc_filter_ordering_prepriorpost.yaml
  testinput/3dvar_qc_filter_ordering_prior.yaml
  testinput/3dvar_qc_filter_ordering_post.yaml
  testinput/3dvar_qc_obserr.yaml
  testinput/3dvar_qc_iterations.yaml
  testinput/3dfgat.yaml
  testinput/4dfgat.yaml
  testinput/4densvar.yaml
  testinput/4densvar_seq.yaml
  testinput/4densvar_hybrid.yaml
  testinput/4densvar_static_crosstime.yaml
  testinput/4dforcing.yaml
  testinput/4dsaddlepoint.yaml
  testinput/4dvar_alpha.yaml
  testinput/4dvar_drgmresr.yaml
  testinput/4dvar_dripcg.yaml
  testinput/4dvar_dripcgqn.yaml
  testinput/4dvar_drpcg.yaml
  testinput/4dvar_drpcgqn.yaml
  testinput/4dvar_drplanczos.yaml
  testinput/4dvar_drplanclmp.yaml
  testinput/4dvar_drplzero.yaml
  testinput/4dvar_fgmres.yaml
  testinput/4dvar_gmresr.yaml
  testinput/4dvar_hybrid.yaml
  testinput/4dvar_ipcg.yaml
  testinput/4dvar_lbgmresr.yaml
  testinput/4dvar_pcg.yaml
  testinput/4dvar_planczos.yaml
  testinput/4dvar_rpcg.yaml
  testinput/4dvar_rplanczos.yaml
  testinput/4dvar_minres.yaml
  testinput/4dvar_modbias.yaml
  testinput/4dvar_obsbias.yaml
  testinput/4dvar_allbiases.yaml
  testinput/addincrement.yaml
  testinput/addincrement_scaled.yaml
  testinput/adjointforecast.yaml
  testinput/diffstates.yaml
  testinput/eda_3dfgat_1.yaml
  testinput/eda_3dfgat_2.yaml
  testinput/eda_3dfgat_3.yaml
  testinput/eda_3dfgat_4.yaml
  testinput/eda_3dfgat.yaml
  testinput/eda_3dvar_1.yaml
  testinput/eda_3dvar_2.yaml
  testinput/eda_3dvar_3.yaml
  testinput/eda_3dvar_4.yaml
  testinput/eda_3dvar.yaml
  testinput/eda_3dvar_block_1.yaml
  testinput/eda_3dvar_block_2.yaml
  testinput/eda_3dvar_block.yaml
  testinput/eda_3dvar_pff_1.yaml
  testinput/eda_3dvar_pff_2.yaml
  testinput/eda_3dvar_pff_3.yaml
  testinput/eda_3dvar_pff_4.yaml
  testinput/eda_3dvar_pff.yaml
  testinput/eda_3dvar_zeromeanpert_1.yaml
  testinput/eda_3dvar_zeromeanpert_2.yaml
  testinput/eda_3dvar_zeromeanpert_3.yaml
  testinput/eda_3dvar_zeromeanpert_4.yaml
  testinput/eda_3dvar_zeromeanpert_compare.yaml
  testinput/eda_3dvar_zeromeanpert.yaml
  testinput/eda_4dvar_3.yaml
  testinput/eda_4dvar.yaml
  testinput/enshofx_1.yaml
  testinput/enshofx_2.yaml
  testinput/enshofx_3.yaml
  testinput/enshofx_4.yaml
  testinput/enshofx.yaml
  testinput/ensvariance.yaml
  testinput/errorcovariance.yaml
  testinput/forecast.yaml
  testinput/forecast_pseudomodel.yaml
  testinput/forecast_identitymodel.yaml
  testinput/fsoi_3dvar_dripcg.yaml
  testinput/fsoi_3dvar_pcg.yaml
  testinput/genenspert.yaml
  testinput/geometry.yaml
  testinput/geometry_iterator.yaml
  testinput/geovals.yaml
  testinput/getkf.yaml
  testinput/getkf_offline_hofx.yaml
  testinput/getkf_offline_hofx_one_member.yaml
  testinput/getvalues.yaml
  testinput/hofx.yaml
  testinput/hofx_tinterp.yaml
  testinput/hofx3d.yaml
  testinput/hofx3d_for_getkf.yaml
  testinput/identitymodel.yaml
  testinput/increment.yaml
  testinput/letkf_gsi.yaml
  testinput/letkf.yaml
  testinput/letkf_noobs.yaml
  testinput/letkf_qc.yaml
  testinput/linearmodel.yaml
  testinput/linearmodelfactory.yaml
  testinput/linobsoperator.yaml
  testinput/localization.yaml
  testinput/makeobs3d.yaml
  testinput/makeobs4d.yaml
  testinput/makeobs4d12h.yaml
  testinput/makeobsbias.yaml
  testinput/makeobspert.yaml
  testinput/model.yaml
  testinput/modelaux.yaml
  testinput/modelauxcontrol.yaml
  testinput/obs_iterator.yaml
  testinput/obsaux.yaml
  testinput/obserrordiagparameters.yaml
  testinput/obslocalization.yaml
  testinput/obslocalizations.yaml
  testinput/obsoperator.yaml
  testinput/obsspace.yaml
  testinput/obsvector.yaml
  testinput/pseudomodel.yaml
  testinput/pseudomodel_state4d.yaml
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
  testinput/sqrtvertloc.yaml
  testinput/state.yaml
  testinput/truth.yaml
)

list( APPEND l95_test_data
  forecast.an.2010-01-01T00:00:00Z.l95
  truth.an.2010-01-01T00:00:00Z.l95
  truth.fc.2010-01-01T00:00:00Z.PT3H.l95
  truth.2010-01-02T00:00:00Z.obt
  truth_static_crosstime.2010-01-01T09:00:00Z.obt
  simplifiedl95.fc.l95
  simplifiedl95.truth.obt
  noobs.obt
  zero.an.2010-01-02T00:00:00Z.l95
)

list( APPEND l95_testoutput
  testoutput/3dvar.test
  testoutput/3dvar_4fsoi_dripcg.test
  testoutput/3dvar_4fsoi_pcg.test
  testoutput/3dvar_noobs.test
  testoutput/3dvar_qc.test
  testoutput/3dvar_qc_filter_ordering_pre.test
  testoutput/3dvar_qc_filter_ordering_prepriorpost.test
  testoutput/3dvar_qc_filter_ordering_prior.test
  testoutput/3dvar_qc_filter_ordering_post.test
  testoutput/3dvar_qc_obserr.test
  testoutput/3dvar_qc_iterations.test
  testoutput/3dfgat.test
  testoutput/4dfgat.test
  testoutput/4densvar.test
  testoutput/4densvar_hybrid.test
  testoutput/4densvar_static_crosstime.test
  testoutput/4dforcing.test
  testoutput/4dsaddlepoint.test
  testoutput/4dvar_alpha.test
  testoutput/4dvar_drgmresr.test
  testoutput/4dvar_dripcg.test
  testoutput/4dvar_dripcgqn.test
  testoutput/4dvar_drpcg.test
  testoutput/4dvar_drpcgqn.test
  testoutput/4dvar_drplanclmp.test
  testoutput/4dvar_drplanczos.test
  testoutput/4dvar_drplzero.test
  testoutput/4dvar_fgmres.test
  testoutput/4dvar_gmresr.test
  testoutput/4dvar_hybrid.test
  testoutput/4dvar_ipcg.test
  testoutput/4dvar_lbgmresr.test
  testoutput/4dvar_modbias.test
  testoutput/4dvar_obsbias.test
  testoutput/4dvar_allbiases.test
  testoutput/4dvar_pcg.test
  testoutput/4dvar_planczos.test
  testoutput/4dvar_rpcg.test
  testoutput/4dvar_rplanczos.test
  testoutput/4dvar_minres.test
  testoutput/diffstates.test
  testoutput/addincrement.test
  testoutput/addincrement_scaled.test
  testoutput/adjointforecast.test
  testoutput/eda_3dfgat.test
  testoutput/eda_3dvar.test
  testoutput/eda_3dvar_block.test
  testoutput/eda_3dvar_pff.test
  testoutput/eda_3dvar_zeromeanpert.test
  testoutput/eda_3dvar_zeromeanpert_compare.test
  testoutput/eda_4dvar.test
  testoutput/enshofx.test
  testoutput/ensvariance.test
  testoutput/forecast.test
  testoutput/forecast_pseudomodel.test
  testoutput/forecast_identitymodel.test
  testoutput/fsoi_3dvar_dripcg.test
  testoutput/fsoi_3dvar_pcg.test
  testoutput/genenspert.test
  testoutput/getkf.test
  testoutput/getkf_offline_hofx.test
  testoutput/getkf_offline_hofx_one_member.test
  testoutput/hofx.test
  testoutput/hofx_tinterp.test
  testoutput/hofx3d.test
  testoutput/hofx3d_for_getkf.test
  testoutput/letkf.test
  testoutput/letkf_noobs.test
  testoutput/letkf_qc.test
  testoutput/letkf_gsi.test
  testoutput/makeobs3d.test
  testoutput/makeobs4d.test
  testoutput/makeobs4d12h.test
  testoutput/makeobsbias.test
  testoutput/makeobspert.test
  testoutput/truth.test
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
  testoutput/sqrtvertloc.test
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

ecbuild_add_test( TARGET test_l95_truth
                  COMMAND  l95_forecast.x
                  ARGS testinput/truth.yaml )

ecbuild_add_test( TARGET test_l95_makeobs4d
                  COMMAND l95_hofx.x
                  ARGS testinput/makeobs4d.yaml
                  TEST_DEPENDS test_l95_truth )


#####################################################################
# interface tests
#####################################################################

ecbuild_add_test( TARGET test_l95_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS "testinput/geometry.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_geometry_iterator
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS "testinput/geometry_iterator.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS "testinput/geovals.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_state
                  SOURCES executables/TestState.cc
                  ARGS "testinput/state.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_getvalues
                  SOURCES executables/TestGetValues.cc
                  ARGS "testinput/getvalues.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_modelauxcontrol
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS "testinput/modelauxcontrol.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_model
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/model.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_pseudomodel
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/pseudomodel.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_pseudomodel_state4d
                  SOURCES executables/TestPseudoModelState4D.cc
                  ARGS "testinput/pseudomodel_state4d.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_identitymodel
                  SOURCES executables/TestModel.cc
                  ARGS "testinput/identitymodel.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS "testinput/increment.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_errorcovariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS "testinput/errorcovariance.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_modelauxincrement
                  SOURCES executables/TestModelAuxIncrement.cc
                  ARGS "testinput/modelaux.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_modelauxcovariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS "testinput/modelaux.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_linearmodel
                  SOURCES executables/TestLinearModel.cc
                  ARGS "testinput/linearmodel.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsspace
                  SOURCES executables/TestObsSpace.cc
                  ARGS "testinput/obsspace.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obs_iterator
                  SOURCES executables/TestObsIterator.cc
                  ARGS "testinput/obs_iterator.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obslocalization
                  SOURCES executables/TestObsLocalization.cc
                  ARGS "testinput/obslocalization.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obslocalizations
                  SOURCES executables/TestObsLocalizations.cc
                  ARGS "testinput/obslocalizations.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsvector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/obsvector.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsdatavector
                  SOURCES executables/TestObsDataVector.cc
                  ARGS "testinput/obsspace.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsoperator
                  SOURCES executables/TestObsOperator.cc
                  ARGS "testinput/obsoperator.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_linobsoperator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS "testinput/linobsoperator.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obserrorcovar
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/obsspace.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsaux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS "testinput/obsaux.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsauxcovariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS "testinput/obsspace.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_obsauxincrement
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS "testinput/obsspace.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS "testinput/localization.yaml"
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_truth )

#####################################################################
# forecast-related tests
#####################################################################

ecbuild_add_test( TARGET test_l95_forecast
                  COMMAND l95_forecast.x
                  ARGS testinput/forecast.yaml )

ecbuild_add_test( TARGET test_l95_forecast_pseudomodel
                  COMMAND l95_forecast.x
                  ARGS testinput/forecast_pseudomodel.yaml )

ecbuild_add_test( TARGET test_l95_forecast_identitymodel
                  COMMAND l95_forecast.x
                  ARGS testinput/forecast_identitymodel.yaml )

#####################################################################
# obs-related tests
#####################################################################

ecbuild_add_test( TARGET test_l95_makeobs3d
                  COMMAND l95_hofx.x
                  ARGS testinput/makeobs3d.yaml
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_makeobsbias
                  COMMAND l95_hofx.x
                  ARGS testinput/makeobsbias.yaml
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_makeobs4d12h
                  COMMAND l95_hofx.x
                  ARGS testinput/makeobs4d12h.yaml
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_makeobspert
                  COMMAND l95_hofx.x
                  ARGS testinput/makeobspert.yaml
                  TEST_DEPENDS test_l95_truth )

ecbuild_add_test( TARGET test_l95_hofx
                  COMMAND l95_hofx.x
                  ARGS testinput/hofx.yaml
                  TEST_DEPENDS test_l95_truth test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_hofx_tinterp
                  COMMAND l95_hofx.x
                  ARGS testinput/hofx_tinterp.yaml
                  TEST_DEPENDS test_l95_truth test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_hofx3d
                  COMMAND l95_hofx3d.x
                  ARGS testinput/hofx3d.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

#####################################################################
# ensemble-related tests
#####################################################################

ecbuild_add_test( TARGET test_l95_genenspert
                  COMMAND l95_genpert.x
                  ARGS testinput/genenspert.yaml )

ecbuild_add_test( TARGET test_l95_enshofx
                  MPI 4
                  COMMAND l95_enshofx.x
                  ARGS testinput/enshofx.yaml
                  DEPENDS l95_enshofx.x
                  TEST_DEPENDS test_l95_genenspert test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_sqrtvertloc
                  COMMAND l95_sqrtofvertloc.x
                  ARGS testinput/sqrtvertloc.yaml )

ecbuild_add_test( TARGET test_l95_ensvariance
                  COMMAND l95_ens_variance.x
                  ARGS testinput/ensvariance.yaml
                  TEST_DEPENDS test_l95_genenspert )

#####################################################################
# 3d variational tests
#####################################################################

ecbuild_add_test( TARGET test_l95_3dvar
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_noobs
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_noobs.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_qc
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_qc.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_qc_filter_ordering_pre
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_qc_filter_ordering_pre.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_qc_filter_ordering_prepriorpost
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_qc_filter_ordering_prepriorpost.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_qc_filter_ordering_prior
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_qc_filter_ordering_prior.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_qc_filter_ordering_post
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_qc_filter_ordering_post.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_qc_obserr
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_qc_obserr.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dvar_qc_iterations
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_qc_iterations.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_3dfgat
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dfgat.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_4dfgat
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dfgat.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )


#####################################################################
# 4d variational tests
#####################################################################

ecbuild_add_test( TARGET test_l95_4dvar_drgmresr
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_drgmresr.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_dripcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_dripcg.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_dripcgqn
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_dripcgqn.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_drpcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_drpcg.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_drpcgqn
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_drpcgqn.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_drplanczos
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_drplanczos.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_drplanclmp
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_drplanclmp.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_fgmres
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_fgmres.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_gmresr
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_gmresr.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_ipcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_ipcg.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_lbgmresr
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_lbgmresr.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_pcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_pcg.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_planczos
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_planczos.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_rpcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_rpcg.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_rplanczos
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_rplanczos.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_minres
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_minres.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_modbias
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_modbias.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_drplzero
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_drplzero.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

#ecbuild_add_test( TARGET test_l95_4dforcing
#                  COMMAND l95_4dvar.x
#                  MPI 2
#                  ARGS testinput/4dforcing.yaml
#                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dsaddlepoint
                  COMMAND l95_4dvar.x
                  MPI 2
                  ARGS testinput/4dsaddlepoint.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_4dvar_obsbias
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_obsbias.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobsbias )

ecbuild_add_test( TARGET test_l95_4dvar_allbiases
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_allbiases.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobsbias )

ecbuild_add_test( TARGET test_l95_4dvar_alpha
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_alpha.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_4dvar_hybrid
                  COMMAND l95_4dvar.x
                  ARGS testinput/4dvar_hybrid.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_4densvar
                  COMMAND l95_4dvar.x
                  MPI 9
                  ARGS testinput/4densvar.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d12h test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_4densvar_seq
                  COMMAND l95_4dvar.x
                  ARGS testinput/4densvar_seq.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d12h test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_4densvar_hybrid
                  COMMAND l95_4dvar.x
                  MPI 9
                  ARGS testinput/4densvar_hybrid.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_4densvar_static_crosstime
                  COMMAND l95_4dvar.x
                  MPI 9
                  ARGS testinput/4densvar_static_crosstime.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d test_l95_genenspert )

#####################################################################
# EDA tests
#####################################################################

ecbuild_add_test( TARGET test_l95_eda_3dfgat
                  MPI 4
                  COMMAND l95_eda.x
                  ARGS testinput/eda_3dfgat.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_eda_3dvar
                  MPI 4
                  COMMAND l95_eda.x
                  ARGS testinput/eda_3dvar.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_eda_4dvar
                  MPI 1
                  COMMAND l95_eda.x
                  ARGS testinput/eda_4dvar.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs4d )

ecbuild_add_test( TARGET test_l95_eda_3dvar_block
                  MPI 2
                  COMMAND l95_eda.x
                  ARGS testinput/eda_3dvar_block.yaml
                  TEST_DEPENDS test_l95_genenspert test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_eda_3dvar_pff
                  MPI 4
                  COMMAND l95_eda.x
                  ARGS testinput/eda_3dvar_pff.yaml
                  TEST_DEPENDS test_l95_genenspert test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_eda_3dvar_zeromeanpert
                  MPI 4
                  COMMAND l95_eda.x
                  ARGS testinput/eda_3dvar_zeromeanpert.yaml
                  DEPENDS l95_eda.x
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_eda_3dvar_zeromeanpert_compare
                  SOURCES lorenz95/ObsErrorDiagZeroMeanPerturbations.cc
                  ARGS testinput/eda_3dvar_zeromeanpert_compare.yaml
                  LIBS lorenz95
                  TEST_DEPENDS test_l95_eda_3dvar_zeromeanpert )
                  
#####################################################################
# FSOI tests
#####################################################################

ecbuild_add_test( TARGET test_l95_3dvar_4fsoi_pcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_4fsoi_pcg.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_fsoi_3dvar_pcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/fsoi_3dvar_pcg.yaml
                  TEST_DEPENDS test_l95_3dvar_4fsoi_pcg )

ecbuild_add_test( TARGET test_l95_3dvar_4fsoi_dripcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/3dvar_4fsoi_dripcg.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_makeobs3d )

ecbuild_add_test( TARGET test_l95_fsoi_3dvar_dripcg
                  COMMAND l95_4dvar.x
                  ARGS testinput/fsoi_3dvar_dripcg.yaml
                  TEST_DEPENDS test_l95_3dvar_4fsoi_dripcg )

#####################################################################
# state-related tests
#####################################################################

ecbuild_add_test( TARGET test_l95_diffstates
                  COMMAND l95_diffstates.x
                  ARGS testinput/diffstates.yaml
                  TEST_DEPENDS test_l95_eda_3dvar test_l95_eda_4dvar )

ecbuild_add_test( TARGET test_l95_addincrement
                  COMMAND l95_addincrement.x
                  ARGS testinput/addincrement.yaml
                  TEST_DEPENDS test_l95_diffstates )

ecbuild_add_test( TARGET test_l95_addincrement_scaled
                  COMMAND l95_addincrement.x
                  ARGS testinput/addincrement_scaled.yaml
                  TEST_DEPENDS test_l95_diffstates )


#####################################################################
# LETKF tests
#####################################################################

ecbuild_add_test( TARGET test_l95_letkf
                  COMMAND l95_letkf.x
                  ARGS testinput/letkf.yaml
                  OMP 2
                  TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_letkf_noobs
                  COMMAND l95_letkf.x
                  ARGS testinput/letkf_noobs.yaml
                  TEST_DEPENDS test_l95_forecast test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_letkf_qc
                  COMMAND l95_letkf.x
                  ARGS testinput/letkf_qc.yaml
                  TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_letkf_gsi
                  COMMAND l95_letkf.x
                  ARGS testinput/letkf_gsi.yaml
                  TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_getkf
                  COMMAND l95_letkf.x
                  ARGS testinput/getkf.yaml
                  TEST_DEPENDS test_l95_makeobs3d test_l95_genenspert )

ecbuild_add_test( TARGET test_l95_hofx3d_for_getkf
                  COMMAND l95_letkf.x
                  ARGS testinput/hofx3d_for_getkf.yaml
                  TEST_DEPENDS test_l95_letkf )

ecbuild_add_test( TARGET test_l95_getkf_offline_hofx
                  COMMAND l95_letkf.x
                  ARGS testinput/getkf_offline_hofx.yaml
                  TEST_DEPENDS test_l95_hofx3d_for_getkf )

ecbuild_add_test( TARGET test_l95_getkf_offline_hofx_one_member
                  COMMAND l95_letkf.x
                  ARGS testinput/getkf_offline_hofx_one_member.yaml
                  TEST_DEPENDS test_l95_hofx3d_for_getkf test_l95_getkf )

#####################################################################
# tests of factory and parameters classes
#####################################################################

ecbuild_add_test( TARGET test_l95_linearmodelfactory
                  SOURCES lorenz95/LinearModelFactory.cc
                  ARGS "testinput/linearmodelfactory.yaml"
                  LIBS lorenz95 )

ecbuild_add_test( TARGET test_l95_obserrordiagparameters
                  SOURCES lorenz95/ObsErrorDiagParameters.cc
                  ARGS "testinput/obserrordiagparameters.yaml"
                  LIBS lorenz95 )

#####################################################################
# tests with simplified L95
#####################################################################

ecbuild_add_test( TARGET test_l95_simplifiedl95_DRGMRESR
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_DRGMRESR.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_DRIPCG
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_DRIPCG.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_DRPCG
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_DRPCG.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_DRPFOM
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_DRPFOM.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_DRPLanczos
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_DRPLanczos.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_FGMRES
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_FGMRES.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_GMRESR
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_GMRESR.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_IPCG
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_IPCG.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_LBGMRESR
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_LBGMRESR.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_MINRES
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_MINRES.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_PCG
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_PCG.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_PLanczos
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_PLanczos.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_RPCG
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_RPCG.yaml )

ecbuild_add_test( TARGET test_l95_simplifiedl95_RPLanczos
                  COMMAND l95_4dvar.x
                  ARGS testinput/simplifiedl95_RPLanczos.yaml )

#####################################################################
# Adjoint Forecast
#####################################################################

ecbuild_add_test( TARGET test_l95_adjointforecast
                  COMMAND l95_adjointforecast.x
                  ARGS testinput/adjointforecast.yaml
                  TEST_DEPENDS test_l95_truth )
