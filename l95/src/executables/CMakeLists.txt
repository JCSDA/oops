ecbuild_add_executable( TARGET  l95_forecast.x
                        SOURCES Forecast.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_adjointforecast.x
                        SOURCES AdjointForecast.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_genpert.x
                        SOURCES GenEnsPertB.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_sqrtofvertloc.x
                        SOURCES SqrtOfVertLoc.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_ens_mean_variance.x
                        SOURCES EnsMeanVariance.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_4dvar.x
                        SOURCES Main4dvar.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_diffstates.x
                        SOURCES DiffStates.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_addincrement.x
                        SOURCES AddIncrement.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_eda.x
                        SOURCES EDA.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_letkf.x
                        SOURCES LETKF.cc
                        LIBS    lorenz95)

ecbuild_add_executable( TARGET  l95_hofx.x
                        SOURCES HofX.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    lorenz95 )

ecbuild_add_executable( TARGET  l95_enshofx.x
                        SOURCES EnsHofX.cc
                        LIBS    lorenz95 )
