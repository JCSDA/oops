find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_subdirectory( src )
add_subdirectory( test )
