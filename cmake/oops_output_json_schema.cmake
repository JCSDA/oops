# (C) Copyright 2020 Met Office UK
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Create a custom target running the application ${APP} with an input YAML file instructing the
# application to output a JSON Schema file describing the structure of input YAML files it accepts.
#
# The JSON Schema file is written to "${CMAKE_BINARY_DIR}/etc/${APP}.schema.json
if(nlohmann_json_FOUND AND nlohmann_json_schema_validator_FOUND)
  function (oops_output_json_schema APP)
    if( oops_SOURCE_DIR )
      set( OOPS_JSON_SCHEMA_CMAKE_DIR "${oops_SOURCE_DIR}/cmake" )
    elseif( oops_CMAKE_DIR )
      set( OOPS_JSON_SCHEMA_CMAKE_DIR "${oops_CMAKE_DIR}" )
    else()
      set( OOPS_JSON_SCHEMA_CMAKE_DIR "${CMAKE_CURRENT_LIST_DIR}" )
    endif()
    set( JSON_SCHEMA_DIR "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_SYSCONFDIR}" )
    file( MAKE_DIRECTORY ${JSON_SCHEMA_DIR} )
    set( DEPS "${APP}" )
    file( GLOB plugin_yml_names "${CMAKE_BINARY_DIR}/share/plugins/*.yml" )
    foreach( plugin_yml_name ${plugin_yml_names} )
      get_filename_component( dep "${plugin_yml_name}" NAME_WLE )
      list( APPEND DEPS "${dep}" )
    endforeach()
    add_custom_target(
      "${APP}.schema.json"
      ALL
      COMMAND "${APP}" "--output-json-schema=${APP}.schema.json"
      DEPENDS "${APP}"
      WORKING_DIRECTORY "${JSON_SCHEMA_DIR}"
    )
    install( FILES "${JSON_SCHEMA_DIR}/${APP}.schema.json" TYPE SYSCONF )
  endfunction()
else()
  function (oops_output_json_schema APP)
    return()
  endfunction()
endif()
