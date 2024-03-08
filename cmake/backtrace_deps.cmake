# (C) Copyright 2021-2024 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# This CMake file tests which stack trace provider libraries are available on
# a particular system. This is needed to set appropriate flags for
# boost stacktrace (used in src/oops/util/signal_trap.cc).
# Thankfully, this eventually may be superseded by C++23's std::stacktrace feature.

# When debugging, unsetting variables triggers re-detection. Ex:
#unset(STACKTRACE_AVAILABLE_none CACHE)
#unset(STACKTRACE_AVAILABLE_libbacktrace CACHE)
#unset(STACKTRACE_AVAILABLE_addr2line CACHE)

find_path(backtrace_header_dir backtrace.h DOC "Path to the backtrace headers")
find_library(backtrace_lib backtrace DOC "Path to the backtrace library")
if ( backtrace_lib )
	cmake_path(GET backtrace_lib PARENT_PATH backtrace_lib_dir)
	cmake_path(GET backtrace_lib EXTENSION LAST_ONLY backtrace_ext)
	set(backtrace_is_static FALSE)
	if (backtrace_ext STREQUAL "a")
		set(backtrace_is_static TRUE)
	endif()
endif()

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

list( APPEND OOPS_STACKTRACE_LIBS_none "")
list( APPEND OOPS_STACKTRACE_DEFS_none -DBOOST_STACKTRACE_USE_NOOP)

list( APPEND OOPS_STACKTRACE_LIBS_default -ldl)
list( APPEND OOPS_STACKTRACE_DEFS_default "")

if( backtrace_lib )
	list( APPEND OOPS_STACKTRACE_LIBS_libbacktrace -ldl -L${backtrace_lib_dir} -lbacktrace)
	list( APPEND OOPS_STACKTRACE_DEFS_libbacktrace -DBOOST_STACKTRACE_USE_BACKTRACE)
	if ( backtrace_header_dir )
		list( APPEND OOPS_STACKTRACE_DEFS_libbacktrace -DBOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE=<${backtrace_header_dir}/backtrace.h>)
	endif()
endif()

list( APPEND OOPS_STACKTRACE_LIBS_addr2line -ldl -lbacktrace)
list( APPEND OOPS_STACKTRACE_DEFS_addr2line -DBOOST_STACKTRACE_USE_ADDR2LINE)

find_program(addr2line_PATH addr2line)
if(addr2line_PATH)
	message( STATUS "Found addr2line at ${addr2line_PATH}." )
	list(APPEND OOPS_STACKTRACE_DEFS_addr2line -DBOOST_STACKTRACE_ADDR2LINE_LOCATION=${addr2line_PATH})
endif()

list( APPEND OOPS_STACKTRACE_POTENTIAL_PROVIDERS addr2line default none )
if( backtrace_lib )
	#list( PREPEND ...) arrived in CMake 3.15. oops depends on CMake 3.12 and above.
	list( INSERT OOPS_STACKTRACE_POTENTIAL_PROVIDERS 0 libbacktrace )
endif()

# Generate the "gnus" and "nognus" configuration variants.
# These reflect variations of the boost definitions that may be required to compile.
foreach ( provider IN LISTS OOPS_STACKTRACE_POTENTIAL_PROVIDERS )
	set( OOPS_STACKTRACE_LIBS_${provider}_nognus ${OOPS_STACKTRACE_LIBS_${provider}})
	set( OOPS_STACKTRACE_DEFS_${provider}_nognus ${OOPS_STACKTRACE_DEFS_${provider}})
	list( APPEND OOPS_STACKTRACE_DEFS_${provider}_nognus -DBOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED)

	set( OOPS_STACKTRACE_LIBS_${provider}_gnus ${OOPS_STACKTRACE_LIBS_${provider}})
	set( OOPS_STACKTRACE_DEFS_${provider}_gnus ${OOPS_STACKTRACE_DEFS_${provider}})
	list( APPEND OOPS_STACKTRACE_DEFS_${provider}_gnus -D_GNU_SOURCE)

	list( APPEND OOPS_STACKTRACE_POTENTIAL_PROVIDERS_2 ${provider} ${provider}_nognus ${provider}_gnus)
endforeach()

# Test each configuration here.
foreach ( provider IN LISTS OOPS_STACKTRACE_POTENTIAL_PROVIDERS_2 )
	list(APPEND CMAKE_REQUIRED_LIBRARIES ${OOPS_STACKTRACE_LIBS_${provider}})
	list(APPEND CMAKE_REQUIRED_DEFINITIONS ${OOPS_STACKTRACE_DEFS_${provider}})
	check_cxx_source_compiles("${stacktracecode}" OOPS_STACKTRACE_AVAILABLE_${provider})
	set(CMAKE_REQUIRED_LIBRARIES ${saved_libraries})
	set(CMAKE_REQUIRED_DEFINITIONS ${saved_defs})
	if ( OOPS_STACKTRACE_AVAILABLE_${provider} )
		list( APPEND OOPS_STACKTRACE_AVAILABLE_PROVIDERS ${provider} )
		break() # Just find the first available provider.
	endif()
endforeach()

set(CMAKE_REQUIRED_LIBRARIES ${saved_libraries})
set(CMAKE_REQUIRED_DEFINITIONS ${saved_defs})
set(CMAKE_REQUIRED_INCLUDES ${saved_incs})

if (NOT OOPS_STACKTRACE_AVAILABLE_PROVIDERS)
	message(FATAL_ERROR "Cannot find a stacktrace provider.")
endif()

list(GET OOPS_STACKTRACE_AVAILABLE_PROVIDERS 0 OOPS_STACKTRACE_PROVIDER)
message( STATUS "Using this provider for stacktraces: ${OOPS_STACKTRACE_PROVIDER}.")

# Small bit of extra logic to ensure that libbacktrace always is linked wherever necessary.
# See src/CMakeLists.txt for usage.
if (OOPS_STACKTRACE_PROVIDER MATCHES "backtrace")
	if (backtrace_is_static)
		add_library(backtrace STATIC IMPORTED)
	else()
		add_library(backtrace SHARED IMPORTED)
	endif()
	set_property(TARGET backtrace PROPERTY
		IMPORTED_LOCATION "${backtrace_lib}")
endif()

