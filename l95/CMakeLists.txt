find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_subdirectory( src )
add_subdirectory( test )

ecbuild_add_test( TARGET l95_coding_norms
                  TYPE SCRIPT
                  COMMAND cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

