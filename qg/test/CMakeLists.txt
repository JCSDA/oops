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
  testinput/ens_variance.yaml
  testinput/forecast.yaml
  testinput/gen_ens_pert_B.yaml
  testinput/hofx.yaml
  testinput/hofx3d.yaml
  testinput/interfaces.yaml
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
  testoutput/ens_variance.test
  testoutput/forecast.test
  testoutput/gen_ens_pert_B.test
  testoutput/hofx.test
  testoutput/hofx3d.test
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

# This line copies all files to binary dir
#file( COPY ${qg_test_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

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

# Create bump directory and links
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_BINARY_DIR}/bump/parametersbump.vbal_vbal_0001-0001.nc
      ${CMAKE_CURRENT_BINARY_DIR}/bump/interfaces_vbal_0001-0001.nc )

#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

ecbuild_add_test( TARGET test_qg_truth
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/truth.yaml"
                       testoutput/truth.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_analytic_forecast
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/analytic_forecast.yaml"
                       testoutput/analytic_forecast.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_forecast
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/forecast.yaml"
                       testoutput/forecast.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_make_obs_3d
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_make_obs.x testinput/make_obs_3d.yaml"
                       testoutput/make_obs_3d.test
                  DEPENDS qg_make_obs.x )

ecbuild_add_test( TARGET test_qg_make_obs_4d_12h
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_make_obs.x testinput/make_obs_4d_12h.yaml"
                       testoutput/make_obs_4d_12h.test
                  DEPENDS qg_make_obs.x )

ecbuild_add_test( TARGET test_qg_make_obs_4d_24h
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_make_obs.x testinput/make_obs_4d_24h.yaml"
                       testoutput/make_obs_4d_24h.test
                  DEPENDS qg_make_obs.x )

ecbuild_add_test( TARGET test_qg_make_obs_4d_biased
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_make_obs.x testinput/make_obs_4d_biased.yaml"
                       testoutput/make_obs_4d_biased.test
                  DEPENDS qg_make_obs.x )

ecbuild_add_test( TARGET test_qg_static_b_init
                  OMP    2
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_staticbinit.x testinput/static_b_init.yaml"
                       testoutput/static_b_init.test
                  DEPENDS qg_staticbinit.x )

# Test QG internal classes

# None yet

# Test interface classes with QG

ecbuild_add_test( TARGET  test_qg_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_geometry_iterator
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_geovals
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_state
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_model_aux
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_model
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_increment
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_error_covariance
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

#ecbuild_add_test( TARGET  test_qg_modelauxincrement
#                  SOURCES executables/TestModelAuxIncrement.cc
#                  ARGS    "testinput/interfaces.yaml"
#                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_model_aux_covariance
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_linear_model
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_locations
                  SOURCES executables/TestLocations.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_observation_space
                  SOURCES executables/TestObservationSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_obs_vector
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obs_operator
                  SOURCES executables/TestObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_lin_obs_operator
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_obs_error_covar
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obs_aux
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obs_aux_increment
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obs_aux_covariance
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_localization
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_linear_variable_change
                  MPI    1
                  SOURCES executables/TestLinearVariableChange.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

# Test oops classes with QG


# Test oops applications with QG

ecbuild_add_test( TARGET test_qg_gen_ens_pert_B
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_gen_ens_pert_B.x testinput/gen_ens_pert_B.yaml"
                       testoutput/gen_ens_pert_B.test
                  DEPENDS qg_gen_ens_pert_B.x )

ecbuild_add_test( TARGET test_qg_hofx
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_hofx.x testinput/hofx.yaml"
                       testoutput/hofx.test
                  DEPENDS qg_hofx.x )

ecbuild_add_test( TARGET test_qg_hofx3d
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_hofx3d.x testinput/hofx3d.yaml"
                       testoutput/hofx3d.test
                  DEPENDS qg_hofx3d.x )

ecbuild_add_test( TARGET test_qg_ens_forecast
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_ens_forecast.x testinput/ens_forecast.yaml"
                       testoutput/ens_forecast.test
                  DEPENDS qg_ens_forecast.x )

ecbuild_add_test( TARGET test_qg_ens_hofx
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_ens_hofx.x testinput/ens_hofx.yaml"
                       testoutput/ens_hofx.test
                  DEPENDS qg_ens_hofx.x )

ecbuild_add_test( TARGET test_qg_ens_variance
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_ens_variance.x testinput/ens_variance.yaml"
                       testoutput/ens_variance.test
                  DEPENDS qg_ens_variance.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_cov
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x testinput/parameters_bump_cov.yaml"
                       testoutput/parameters_bump_cov.test
                  DEPENDS qg_estimate_parameters.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_lct
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x testinput/parameters_bump_lct.yaml"
                       testoutput/parameters_bump_lct.test
                  DEPENDS qg_estimate_parameters.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_hyb
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x testinput/parameters_bump_hyb.yaml"
                       testoutput/parameters_bump_hyb.test
                  DEPENDS qg_estimate_parameters.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_loc_3d
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x testinput/parameters_bump_loc_3d.yaml"
                       testoutput/parameters_bump_loc_3d.test
                  DEPENDS qg_estimate_parameters.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_loc_4d
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_estimate_parameters.x testinput/parameters_bump_loc_4d.yaml"
                       testoutput/parameters_bump_loc_4d.test
                  DEPENDS qg_estimate_parameters.x )

ecbuild_add_test( TARGET test_qg_bump_links
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/bump_links.sh )

ecbuild_add_test( TARGET test_qg_dirac_bump_cov
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_bump_cov.yaml"
                       testoutput/dirac_bump_cov.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_bump_hyb
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_bump_hyb.yaml"
                       testoutput/dirac_bump_hyb.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_bump_lct
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_bump_lct.yaml"
                       testoutput/dirac_bump_lct.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_bump_loc_3d
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_bump_loc_3d.yaml"
                       testoutput/dirac_bump_loc_3d.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_bump_loc_4d
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_bump_loc_4d.yaml"
                       testoutput/dirac_bump_loc_4d.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_qg_cov
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_qg_cov.yaml"
                       testoutput/dirac_qg_cov.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_qg_hyb
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_qg_hyb.yaml"
                       testoutput/dirac_qg_hyb.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_qg_loc_3d
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_qg_loc_3d.yaml"
                       testoutput/dirac_qg_loc_3d.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_qg_loc_4d
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_qg_loc_4d.yaml"
                       testoutput/dirac_qg_loc_4d.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_no_loc
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac_no_loc.yaml"
                       testoutput/dirac_no_loc.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_3densvar
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3densvar.yaml"
                       testoutput/3densvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3densvar_bump
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3densvar_bump.yaml"
                       testoutput/3densvar_bump.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dvar
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar.yaml"
                       testoutput/3dvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dvar_bump
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar_bump.yaml"
                       testoutput/3dvar_bump.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dvar_change_var
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar_change_var.yaml"
                       testoutput/3dvar_change_var.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dvar_hybrid
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar_hybrid.yaml"
                       testoutput/3dvar_hybrid.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dvar_hybrid_bump
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar_hybrid_bump.yaml"
                       testoutput/3dvar_hybrid_bump.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dfgat
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dfgat.yaml"
                       testoutput/3dfgat.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4densvar
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4densvar.yaml"
                       testoutput/4densvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4densvar_bump
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4densvar_bump.yaml"
                       testoutput/4densvar_bump.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_dripcg
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_dripcg.yaml"
                       testoutput/4dvar_dripcg.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drpcg_lmp
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_drpcg_lmp.yaml"
                       testoutput/4dvar_drpcg_lmp.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drpfom
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_drpfom.yaml"
                       testoutput/4dvar_drpfom.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_drplanczos.yaml"
                       testoutput/4dvar_drplanczos.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_bump
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_drplanczos_bump.yaml"
                       testoutput/4dvar_drplanczos_bump.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_drplanczos_hybrid.yaml"
                       testoutput/4dvar_drplanczos_hybrid.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos_hybrid_bump
                  MPI    1
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_drplanczos_hybrid_bump.yaml"
                       testoutput/4dvar_drplanczos_hybrid_bump.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_forcing
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_forcing.yaml"
                       testoutput/4dvar_forcing.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_ipcg
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_ipcg.yaml"
                       testoutput/4dvar_ipcg.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_obs_biased
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_obs_biased.yaml"
                       testoutput/4dvar_obs_biased.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_rpcg
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_rpcg.yaml"
                       testoutput/4dvar_rpcg.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_saddlepoint
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar_saddlepoint.yaml"
                       testoutput/4dvar_saddlepoint.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_convertstate
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_convertstate.x
                  ARGS testinput/convertstate.yaml
                  DEPENDS qg_convertstate.x )

ecbuild_add_test( TARGET test_qg_dfi
                  OMP    2
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dfi.x testinput/dfi.yaml"
                       testoutput/dfi.test
                  DEPENDS qg_dfi.x )

ecbuild_add_test( TARGET test_qg_eda_4dvar_pert_obs
                  MPI    4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/eda_4dvar_pert_obs.yaml
                  DEPENDS qg_eda.x )

#ecbuild_add_test( TARGET test_qg_eda_4dvar_pert_obs
#                  MPI    4
#                  TYPE SCRIPT
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
#                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_eda.x testinput/eda_4dvar_pert_obs.yaml"
#                       testoutput/eda_4dvar_pert_obs.test
#                  DEPENDS qg_eda.x )
