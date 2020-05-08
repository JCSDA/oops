function(oops_add_test)
  set( options COMPARE )
  set( single_value_args TESTNAME MODELNAME YAMLNAME EXENAME MPI OMP CTOL IDIF RUN_FILE REF_FILE )
  set( multi_value_args  DEPENDS TEST_DEPENDS)
  cmake_parse_arguments( _p "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  # Set default values

  if (_p_MPI)
    if( NOT MPIEXEC )
      if( MPIEXEC_EXECUTABLE )
        set( MPIEXEC ${MPIEXEC_EXECUTABLE} )
      else()
        find_program( MPIEXEC NAMES mpiexec mpirun lamexec srun
                    DOC "Executable for running MPI programs." )
      endif()
    endif()
    set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC")
    set(MPI_cmd ${MPIEXEC}\ ${MPIEXEC_NUMPROC_FLAG}\ ${_p_MPI})
  else()
    set (MPI_cmd "")
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
                         ${MPI_cmd}
                    OMP ${_p_OMP}
                    DEPENDS ${_p_DEPENDS}
                    TEST_DEPENDS ${_p_TEST_DEPENDS}
                   )
endfunction(oops_add_test)
