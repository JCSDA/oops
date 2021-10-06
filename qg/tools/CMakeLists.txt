# Link scripts
list( APPEND test_files
    qg_cost.py
    qg_fields.py
    qg_obs.py
)

foreach( FILENAME ${test_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                             ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
                             ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach()
