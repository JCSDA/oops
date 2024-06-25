list( APPEND suppression_files
    lsan.supp
    valgrind.supp
)

set( destination share/oops/suppressions )

install(
    FILES       ${suppression_files}
    DESTINATION ${destination}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${destination} )

foreach( file ${suppression_files} )
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "copy_if_different"
        "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
        "${CMAKE_BINARY_DIR}/${destination}/${file}")
endforeach()
