# (C) Copyright 2020 Met Office UK
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set(OOPS_JSON_SCHEMA_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})

# Create a custom target running the application ${APP} with an input YAML file instructing the
# application to output a JSON Schema file describing the structure of input YAML files it accepts.
#
# The JSON Schema file is written to "${CMAKE_BINARY_DIR}/etc/${APP}.schema.json
function (oops_output_json_schema APP)
  set( JSON_SCHEMA_DIR "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_SYSCONFDIR}" )
  set( JSON_SCHEMA_GENERATOR_INPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/schemageneratorinputs" )
  file( MAKE_DIRECTORY ${JSON_SCHEMA_DIR} )
  file( MAKE_DIRECTORY ${JSON_SCHEMA_GENERATOR_INPUT_DIR} )

  set( JSON_SCHEMA_PATH "${JSON_SCHEMA_DIR}/${APP}.schema.json" ) # used by configure_file below
  configure_file( "${OOPS_JSON_SCHEMA_CMAKE_DIR}/output_json_schema.yaml.in"
                  "${JSON_SCHEMA_GENERATOR_INPUT_DIR}/${APP}.yaml" )
  configure_file( "${OOPS_JSON_SCHEMA_CMAKE_DIR}/oops_output_json_schema_x.cmake.in"
                  "${JSON_SCHEMA_GENERATOR_INPUT_DIR}/${APP}.cmake" )
  add_custom_command( OUTPUT "${JSON_SCHEMA_DIR}/${APP}.schema.json"
                      DEPENDS "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${APP}"
                              "${JSON_SCHEMA_GENERATOR_INPUT_DIR}/${APP}.yaml"
                              "${JSON_SCHEMA_GENERATOR_INPUT_DIR}/${APP}.cmake"
                      COMMAND "${CMAKE_COMMAND}"
                              -P "${JSON_SCHEMA_GENERATOR_INPUT_DIR}/${APP}.cmake" )
  add_custom_target( "generate_${APP}.schema.json"
                     ALL
                     DEPENDS "${JSON_SCHEMA_DIR}/${APP}.schema.json" )
endfunction()
