# (C) Copyright 2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# This CMake file tests which stack trace provider libraries are available on
# a particular system. This is needed to set appropriate flags for
# boost stacktrace (used in src/oops/util/signal_trap.cc).

include(CheckCXXSourceCompiles)

# check_cxx_source_compiles uses global flags, unfortunately, so we need to
# save a bit of state to run the tests.

set(saved_libraries ${CMAKE_REQUIRED_LIBRARIES})
set(saved_defs ${CMAKE_REQUIRED_DEFINITIONS})
set(saved_incs ${CMAKE_REQUIRED_INCLUDES})

list( APPEND CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIRS})

string(CONFIGURE [[
    #include <iostream>
    #include <boost/stacktrace.hpp>

    int main() {
        std::cout << boost::stacktrace::stacktrace() << std::endl;
        return 0;
    }
]] stacktracecode @ONLY)


# Different configs

list( APPEND OOPS_STACKTRACE_none_LIBS "")
list( APPEND OOPS_STACKTRACE_none_DEFS -DBOOST_STACKTRACE_USE_NOOP)

list( APPEND OOPS_STACKTRACE_default_LIBS -ldl)
list( APPEND OOPS_STACKTRACE_default_DEFS "")

list( APPEND OOPS_STACKTRACE_libbacktrace_LIBS -ldl -lbacktrace)
list( APPEND OOPS_STACKTRACE_libbacktrace_DEFS -DBOOST_STACKTRACE_USE_BACKTRACE)

list( APPEND OOPS_STACKTRACE_addr2line_LIBS -ldl -lbacktrace)
list( APPEND OOPS_STACKTRACE_addr2line_DEFS -DBOOST_STACKTRACE_USE_ADDR2LINE)
find_program(addr2line_PATH addr2line)
if(addr2line_PATH)
	message( STATUS "Found addr2line at ${addr2line_PATH}." )
	list(APPEND OOPS_STACKTRACE_addr2line_DEFS -DBOOST_STACKTRACE_ADDR2LINE_LOCATION=${addr2line_PATH})
endif()


# Test each configuration here.
foreach ( provider IN ITEMS libbacktrace addr2line default none )
	list(APPEND CMAKE_REQUIRED_LIBRARIES ${OOPS_STACKTRACE_${provider}_LIBS})
	list(APPEND CMAKE_REQUIRED_DEFINITIONS ${OOPS_STACKTRACE_${provider}_DEFS})
	check_cxx_source_compiles("${stacktracecode}" OOPS_STACKTRACE_${provider}_AVAILABLE)
	set(CMAKE_REQUIRED_LIBRARIES ${saved_libraries})
	set(CMAKE_REQUIRED_DEFINITIONS ${saved_defs})
	if ( OOPS_STACKTRACE_${provider}_AVAILABLE )
		list( APPEND OOPS_STACKTRACE_AVAILABLE_PROVIDERS ${provider} )
	endif()
endforeach()

set(CMAKE_REQUIRED_LIBRARIES ${saved_libraries})
set(CMAKE_REQUIRED_DEFINITIONS ${saved_defs})
set(CMAKE_REQUIRED_INCLUDES ${saved_incs})

if (NOT OOPS_STACKTRACE_AVAILABLE_PROVIDERS)
	message(FATAL_ERROR "Cannot find a stacktrace provider.")
endif()

message( STATUS "Boost stacktrace supports these providers: ${OOPS_STACKTRACE_AVAILABLE_PROVIDERS}.")
list(GET OOPS_STACKTRACE_AVAILABLE_PROVIDERS 0 OOPS_STACKTRACE_PROVIDER)
message( STATUS "Using this provider for stacktraces: ${OOPS_STACKTRACE_PROVIDER}.")


