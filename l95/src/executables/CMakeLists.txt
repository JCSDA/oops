ecbuild_add_executable( TARGET  l95_forecast.x
                        SOURCES Forecast.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_genpert.x
                        SOURCES GenEnsPertB.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_4dvar.x
                        SOURCES Main4dvar.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_makeobs.x
                        SOURCES MakeObs.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_hofx.x
                        SOURCES HofX.cc
                        LIBS    lorenz95 )
