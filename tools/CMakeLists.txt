
# Create Data directory for test input data
list( APPEND tool_files
    compare.py
    cpplint.py
    plot.py
    plot/cost.py
    plot/l95_cost.py
    plot/l95_fields.py
    plot/qg_cost.py
    plot/qg_fields.py
    plot/qg_obs.py
    test_wrapper.sh
)

foreach(FILENAME IN LISTS tool_files)
    set(_tool ${PROJECT_NAME}_${FILENAME})
    # Some entries live in a subdirectory (plot/*), so both copies below have to
    # account for it. `cmake -E copy` used to create the destination directory
    # implicitly; file(COPY_FILE) does not.
    get_filename_component(_tool_subdir ${FILENAME} DIRECTORY)

    #Install in build-tree
    if( _tool_subdir )
      file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_${_tool_subdir} )
    endif()
    file( COPY_FILE
        ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
        ${CMAKE_BINARY_DIR}/bin/${_tool} ONLY_IF_DIFFERENT )
    #Install in install-tree
    install(PROGRAMS ${FILENAME} RENAME ${_tool} DESTINATION bin )

    #TODO remove this once no packages depend on these tools in their original location.
    if( _tool_subdir )
      file( MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/${_tool_subdir} )
    endif()
    file( COPY_FILE
        ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
        ${CMAKE_BINARY_DIR}/bin/${FILENAME} ONLY_IF_DIFFERENT )
endforeach()

