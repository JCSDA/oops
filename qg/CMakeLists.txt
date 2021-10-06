
add_subdirectory( model )
add_subdirectory( mains )
add_subdirectory( test )
add_subdirectory( tools )

ecbuild_add_test( TARGET qg_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/oops_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

