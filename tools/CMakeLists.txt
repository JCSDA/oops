
# Create Data directory for test input data
list( APPEND oops_tool_files
  compare.sh
  compare.py
  cpplint.py
)

foreach(FILENAME ${oops_tool_files})
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
     ${PROJECT_SOURCE_DIR}/tools/${FILENAME}
     ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach(FILENAME)

install(PROGRAMS ${oops_tool_files} DESTINATION bin)
