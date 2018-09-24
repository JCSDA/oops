# first we can indicate the documentation build as an option and
# set it to OFF by default since documenation will usually be
# built at the bundle level
option(BUILD_OOPS_DOC "Build OOPS documentation" OFF)

if (BUILD_OOPS_DOC)

    message("OOPS Documentation Requested")
    
    # check if Doxygen is installed
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen_oops ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating OOPS documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
        message("Doxygen needs to be installed to generate the OOPS documentation")
    endif (DOXYGEN_FOUND)
endif (BUILD_OOPS_DOC)
