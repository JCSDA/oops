#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

ecbuild_add_test( TARGET  test_util_datetime
                  SOURCES util/DateTime.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_duration
                  SOURCES util/Duration.cc
                  LIBS    oops eckit )

if( HAVE_FCTEST )
add_fctest( TARGET  test_util_config_intfb.x
            SOURCES util/config.F90
            CONDITION HAVE_FCTEST
            LIBS    oops )

add_fctest( TARGET  test_util_datetime_intfb.x
            SOURCES util/datetime.F90
            CONDITION HAVE_FCTEST
            LIBS    oops )

add_fctest( TARGET  test_util_duration_intfb.x
            SOURCES util/duration.F90
            CONDITION HAVE_FCTEST
            LIBS    oops )
endif()

