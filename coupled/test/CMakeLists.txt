list( APPEND coupled_testinput
  testinput/coupledconvertstate.yaml
  testinput/coupledgeometry.yaml
  testinput/coupledgeometryparallel.yaml
  testinput/coupledmodel.yaml
  testinput/coupledmodelauxcontrol.yaml
  testinput/coupledmodelparallel.yaml
  testinput/coupledstate.yaml
  testinput/coupledstateparallel.yaml
  testinput/forecast_qg_l95.yaml
)

list( APPEND coupled_testref
  testref/coupledconvertstate.test
  testref/coupledgeometry.test
  testref/coupledmodel.test
  testref/coupledstate.test
  testref/forecast_qg_l95.test
)

# Create directory for test yamls and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${coupled_testinput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create directory for test references and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
foreach(FILENAME ${coupled_testref})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for test output data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

ecbuild_add_resources( TARGET coupled_test_scripts
                       SOURCES_PACK
                       ${coupled_testinput}
                     )

# Test CoupledGeometry class
ecbuild_add_test( TARGET  test_coupled_geometry_qg_l95
                  SOURCES executables/TestCoupledGeometry.cc
                  ARGS    "testinput/coupledgeometry.yaml"
                  LIBS    qg lorenz95 )

# Test CoupledGeometry class in parallel
ecbuild_add_test( TARGET  test_coupled_geometry_qg_l95_parallel
                  SOURCES executables/TestCoupledGeometry.cc
                  ARGS    "testinput/coupledgeometryparallel.yaml"
                  MPI     2
                  LIBS    qg lorenz95 )

# Test CoupledModelAuxControl class
ecbuild_add_test( TARGET  test_coupled_modelauxcontrol_qg_l95
                  SOURCES executables/TestCoupledModelAuxControl.cc
                  ARGS    "testinput/coupledmodelauxcontrol.yaml"
                  LIBS    qg lorenz95 )

# Test CoupledState class
ecbuild_add_test( TARGET  test_coupled_state_qg_l95
                  SOURCES executables/TestCoupledState.cc
                  ARGS    "testinput/coupledstate.yaml"
                  LIBS    qg lorenz95
                  TEST_DEPENDS test_qg_truth test_l95_3dvar)

# Test CoupledState class in parallel
ecbuild_add_test( TARGET  test_coupled_state_qg_l95_parallel
                  SOURCES executables/TestCoupledState.cc
                  ARGS    "testinput/coupledstateparallel.yaml"
                  MPI     2
                  LIBS    qg lorenz95
                  TEST_DEPENDS test_qg_truth test_l95_3dvar)

# Test CoupledVariableChange class
ecbuild_add_test( TARGET  test_coupled_convertstate_qg_l95
                  SOURCES executables/TestCoupledConvertState.cc
                  ARGS    "testinput/coupledconvertstate.yaml"
                  LIBS    qg lorenz95
                  TEST_DEPENDS test_qg_truth test_l95_3dvar)

# Test CoupledModel class
ecbuild_add_test( TARGET  test_coupled_model_qg_l95
                  SOURCES executables/TestCoupledModel.cc
                  ARGS    "testinput/coupledmodel.yaml"
                  LIBS    qg lorenz95
                  TEST_DEPENDS test_qg_truth test_l95_3dvar)

# Test CoupledModel class
ecbuild_add_test( TARGET  test_coupled_model_qg_l95_parallel
                  SOURCES executables/TestCoupledModel.cc
                  ARGS    "testinput/coupledmodelparallel.yaml"
                  MPI     2
                  LIBS    qg lorenz95
                  TEST_DEPENDS test_qg_truth test_l95_3dvar)

# Test CoupledForecast application
ecbuild_add_test( TARGET  test_coupled_forecast_qg_l95
                  SOURCES executables/forecast_qg_l95.cc
                  ARGS    "testinput/forecast_qg_l95.yaml"
                  LIBS    qg lorenz95
                  TEST_DEPENDS test_qg_truth test_l95_3dvar)
