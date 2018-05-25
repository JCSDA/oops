include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/model )

add_subdirectory( model )
add_subdirectory( mains )
add_subdirectory( test )

ecbuild_add_test( TARGET qg_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

