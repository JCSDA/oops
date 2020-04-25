function(oops_add_test)
  set( options COMPARE )
  set( single_value_args TESTNAME MODELNAME YAMLNAME EXENAME MPI OMP CTOL IDIF RUN_FILE REF_FILE )
  set( multi_value_args  DEPENDS TEST_DEPENDS)
  cmake_parse_arguments( _p "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )
 
  # Set default values

  if (NOT _p_MPI)
    set(_p_MPI2 0)
  else() 
    set(_p_MPI2 _p_MPI)
  endif()

  if (NOT _p_RUN_FILE)
    set( _p_RUN_FILE ${_p_TESTNAME}.log.out)
  endif()

  if (NOT _p_REF_FILE)
    set( _p_REF_FILE ${_p_TESTNAME}.test)
  endif()

  if (NOT _p_CTOL)
    set( _p_CTOL 0.0)  # Max relative diff
  endif()

  if (NOT _p_IDIF)
    set( _p_IDIF 0)  # Max diff
  endif()

  ecbuild_add_test( TARGET test_${_p_MODELNAME}_${_p_TESTNAME}
                    TYPE SCRIPT
                    COMMAND ${CMAKE_BINARY_DIR}/bin/test_wrapper.sh
                    ARGS ${CMAKE_BINARY_DIR}/bin/${_p_EXENAME} 
                         ${_p_YAMLNAME} 
                         ${CMAKE_BINARY_DIR}/bin/compare.py 
                         ${_p_RUN_FILE} 
                         ${_p_REF_FILE}
                         ${_p_CTOL}
                         ${_p_IDIF}
                         ${_p_MPI2}
                    OMP ${_p_OMP}
                    MPI ${_p_MPI} 
                    DEPENDS ${_p_DEPENDS} 
                    TEST_DEPENDS ${_p_TEST_DEPENDS} 
                   )
endfunction(oops_add_test)
