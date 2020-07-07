add_subdirectory( src/lorenz95 )
add_subdirectory( src/executables )
add_subdirectory( test )

ecbuild_add_test( TARGET l95_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/oops_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

