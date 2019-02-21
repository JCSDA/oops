list( APPEND qg_test_input
  testinput/3dvar.yaml
  testinput/3dvar.chvar.yaml
  testinput/3dfgat.yaml
  testinput/4densvar.yaml
  testinput/4dvar.alpha.yaml
  testinput/4dvar.bumpcov.yaml
  testinput/4dforcing.yaml
  testinput/4dhybrid.yaml
  testinput/4dhybridbump.yaml
  testinput/4dsaddlepoint.yaml
  testinput/4dvar.dripcg.yaml
  testinput/4dvar.drpcglmp.yaml
  testinput/4dvar.drpfom.yaml
  testinput/4dvar.drplanczos.yaml
  testinput/4dvar.eda.pertobs.yaml
  testinput/4dvar.ipcg.yaml
  testinput/4dvar.obsbias.yaml
  testinput/4dvar.rpcg.yaml
  testinput/dfi.yaml
  testinput/dirac.bumpcov.yaml
  testinput/dirac.qgcov.yaml
  testinput/dirac.bumploc.yaml
  testinput/enshofx.yaml
  testinput/parametersbump.cov.yaml
  testinput/parametersbump.lct.yaml
  testinput/parametersbump.loc.yaml
  testinput/parametersbump.readloc.yaml
  testinput/parametersbump.vbal.yaml
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
  testinput/analytic_init_fc.yaml
  testinput/static_binit.yml
  testinput/convertstate.yaml
)

list( APPEND qg_testoutput
  testoutput/3dvar.test
  testoutput/3dfgat.test
  testoutput/4densvar.test
  testoutput/4dvar.alpha.test
  testoutput/4dvar.bumpcov.test
  testoutput/4dforcing.test
  testoutput/4dhybrid.test
  testoutput/4dhybridbump.test
  testoutput/4dsaddlepoint.test
  testoutput/4dvar.dripcg.test
  testoutput/4dvar.drpcglmp.test
  testoutput/4dvar.drpfom.test
  testoutput/4dvar.drplanczos.test
  testoutput/4dvar.ipcg.test
  testoutput/4dvar.obsbias.test
  testoutput/4dvar.rpcg.test
  testoutput/dfi.test
  testoutput/dirac.bumpcov.test
  testoutput/dirac.qgcov.test
  testoutput/dirac.bumploc.test
  testoutput/enshofx.test
  testoutput/parametersbump.cov.test
  testoutput/parametersbump.lct.test
  testoutput/parametersbump.loc.test
  testoutput/parametersbump.readloc.test
  testoutput/parametersbump.vbal.test
  testoutput/forecast.test
  testoutput/genenspert.test
  testoutput/hofx.test
  testoutput/makeobs3d.test
  testoutput/makeobs4d.test
  testoutput/makeobs4d12h.test
  testoutput/makeobsbias.test
  testoutput/truth.test
  testoutput/analytic_init_fc.test
)

# This line copies all files to binary dir
#file( COPY ${qg_test_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${qg_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for reference output and symlink all files
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

#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

ecbuild_add_test( TARGET test_qg_truth
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/truth.yaml"
                       testoutput/truth.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_aninit_fc
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/analytic_init_fc.yaml"
                       testoutput/analytic_init_fc.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_forecast
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_forecast.x testinput/forecast.yaml"
                       testoutput/forecast.test
                  DEPENDS qg_forecast.x )

ecbuild_add_test( TARGET test_qg_makeobs3d
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobs3d.yaml"
                       testoutput/makeobs3d.test
                  DEPENDS qg_makeobs.x )

ecbuild_add_test( TARGET test_qg_makeobs4d
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobs4d.yaml"
                       testoutput/makeobs4d.test
                  DEPENDS qg_makeobs.x )

ecbuild_add_test( TARGET test_qg_makeobs4d12h
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobs4d12h.yaml"
                       testoutput/makeobs4d12h.test
                  DEPENDS qg_makeobs.x )

ecbuild_add_test( TARGET test_qg_makeobsbias
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_makeobs.x testinput/makeobsbias.yaml"
                       testoutput/makeobsbias.test
                  DEPENDS qg_makeobs.x )

ecbuild_add_test( TARGET  test_qg_staticbinit
                  TYPE    EXE
                  SOURCES ../mains/qgStaticBInit.cc
                  ARGS    "testinput/static_binit.yml"
                  MPI     1		  
                  LIBS    qg )

# Test QG internal classes

# None yet

# Test interface classes with QG

#ecbuild_add_test( TARGET  test_qg_geometry
#                  BOOST
#                  SOURCES executables/TestGeometry.cc
#                  ARGS    "testinput/interfaces.yaml"
#                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_geometry
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.yaml"
                  DEFINITIONS OOPS_TEST_NO_BOOST
                  LIBS    qg )


ecbuild_add_test( TARGET  test_qg_geometry_iterator
                  BOOST
                  SOURCES executables/TestGeometryIterator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_geovals
                  BOOST
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_modelaux
                  BOOST
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

#ecbuild_add_test( TARGET  test_qg_modelauxincrement
#                  BOOST
#                  SOURCES executables/TestModelAuxIncrement.cc
#                  ARGS    "testinput/interfaces.yaml"
#                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_modelauxcovariance
                  BOOST
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_linearmodel
                  BOOST
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_locations
                  BOOST
                  SOURCES executables/TestLocations.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_observationspace
                  BOOST
                  SOURCES executables/TestObservationSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_obsvector
                  BOOST
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsoperator
                  BOOST
                  SOURCES executables/TestObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_linobsoperator
                  BOOST
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET test_qg_obserrorcovar
                  BOOST
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsaux
                  BOOST
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsauxincrement
                  BOOST
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_obsauxcovariance
                  BOOST
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

ecbuild_add_test( TARGET  test_qg_localization
                  BOOST
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

# Application test_qg_linearvariablechange requires test_qg_parameters_bump_vbal output, which requires test_qg_genenspert output
# Vertical balance written by test_qg_parameters_bump_vbal is linked by test_qg_copy_vbal

ecbuild_add_test( TARGET test_qg_genenspert
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_genpert.x testinput/genenspert.yaml"
                       testoutput/genenspert.test
                  DEPENDS qg_genpert.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_vbal
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_parameters.x testinput/parametersbump.vbal.yaml"
                       testoutput/parametersbump.vbal.test
                  DEPENDS qg_parameters.x )

ecbuild_add_test( TARGET test_qg_copy_vbal
                  TYPE SCRIPT
                  COMMAND "ln"
                  ARGS  -sf bump_parametersbump.vbal_vbal_0001-0001.nc bump_interfaces_vbal_0001-0001.nc )

ecbuild_add_test( TARGET  test_qg_linearvariablechange
                  MPI    1
                  BOOST
                  SOURCES executables/TestLinearVariableChange.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    qg )

# Test oops classes with QG


# Test oops applications with QG

ecbuild_add_test( TARGET test_qg_hofx
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_hofx.x testinput/hofx.yaml"
                       testoutput/hofx.test
                  DEPENDS qg_hofx.x )

ecbuild_add_test( TARGET test_qg_enshofx
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_enshofx.x testinput/enshofx.yaml"
                       testoutput/enshofx.test
                  DEPENDS qg_enshofx.x )

ecbuild_add_test( TARGET test_qg_3dvar
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar.yaml"
                       testoutput/3dvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dvar_chvar
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dvar.chvar.yaml"
                       testoutput/3dvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_3dfgat
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/3dfgat.yaml"
                       testoutput/3dfgat.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_bump_covariance
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.bumpcov.yaml"
                       testoutput/4dvar.bumpcov.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drpcglmp
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.drpcglmp.yaml"
                       testoutput/4dvar.drpcglmp.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drpfom
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.drpfom.yaml"
                       testoutput/4dvar.drpfom.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_drplanczos
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.drplanczos.yaml"
                       testoutput/4dvar.drplanczos.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_dripcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.dripcg.yaml"
                       testoutput/4dvar.dripcg.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_ipcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.ipcg.yaml"
                       testoutput/4dvar.ipcg.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_rpcg
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.rpcg.yaml"
                       testoutput/4dvar.rpcg.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_obsbias
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.obsbias.yaml"
                       testoutput/4dvar.obsbias.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_forcing
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dforcing.yaml"
                       testoutput/4dforcing.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_saddlepoint
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dsaddlepoint.yaml"
                       testoutput/4dsaddlepoint.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_hybrid
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dhybrid.yaml"
                       testoutput/4dhybrid.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4dvar_hybrid_bump
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dhybridbump.yaml"
                       testoutput/4dhybridbump.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_4densvar
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4densvar.yaml"
                       testoutput/4densvar.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_alpha
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_4dvar.x testinput/4dvar.alpha.yaml"
                       testoutput/4dvar.alpha.test
                  DEPENDS qg_4dvar.x )

ecbuild_add_test( TARGET test_qg_dirac_bumpcov
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac.bumpcov.yaml"
                       testoutput/dirac.bumpcov.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_bumploc
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac.bumploc.yaml"
                       testoutput/dirac.bumploc.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_dirac_qgcov
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dirac.x testinput/dirac.qgcov.yaml"
                       testoutput/dirac.qgcov.test
                  DEPENDS qg_dirac.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_cov
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_parameters.x testinput/parametersbump.cov.yaml"
                       testoutput/parametersbump.cov.test
                  DEPENDS qg_parameters.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_lct
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_parameters.x testinput/parametersbump.lct.yaml"
                       testoutput/parametersbump.lct.test
                  DEPENDS qg_parameters.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_loc
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_parameters.x testinput/parametersbump.loc.yaml"
                       testoutput/parametersbump.loc.test
                  DEPENDS qg_parameters.x )

ecbuild_add_test( TARGET test_qg_parameters_bump_readloc
                  MPI    1
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_parameters.x testinput/parametersbump.readloc.yaml"
                       testoutput/parametersbump.readloc.test
                  DEPENDS qg_parameters.x )

ecbuild_add_test( TARGET test_qg_dfi
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS "${CMAKE_BINARY_DIR}/bin/qg_dfi.x testinput/dfi.yaml"
                       testoutput/dfi.test
                  DEPENDS qg_dfi.x )

ecbuild_add_test( TARGET test_qg_4dvar_eda
                  MPI 4
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_eda.x
                  ARGS testinput/4dvar.eda.pertobs.yaml
                  DEPENDS qg_eda.x )

ecbuild_add_test( TARGET test_qg_convertstate
                  COMMAND ${CMAKE_BINARY_DIR}/bin/qg_convertstate.x
                  ARGS testinput/convertstate.yaml
                  DEPENDS qg_convertstate.x )
