include(oops_output_json_schema)

ecbuild_add_executable( TARGET  l95_forecast.x
                        SOURCES Forecast.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_forecast.x )

ecbuild_add_executable( TARGET  l95_adjointforecast.x
                        SOURCES AdjointForecast.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_adjointforecast.x )

ecbuild_add_executable( TARGET  l95_genpert.x
                        SOURCES GenEnsPertB.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_genpert.x )

ecbuild_add_executable( TARGET  l95_sqrtofvertloc.x
                        SOURCES SqrtOfVertLoc.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_sqrtofvertloc.x )

ecbuild_add_executable( TARGET  l95_ens_variance.x
                        SOURCES EnsVariance.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_ens_variance.x )

ecbuild_add_executable( TARGET  l95_ens_mean_variance.x
                        SOURCES EnsMeanVariance.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_ens_mean_variance.x )

ecbuild_add_executable( TARGET  l95_4dvar.x
                        SOURCES Main4dvar.cc
                        LIBS    lorenz95 )
# TODO: oops_output_json_schema( l95_4dvar.x )

ecbuild_add_executable( TARGET  l95_diffstates.x
                        SOURCES DiffStates.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_diffstates.x )

ecbuild_add_executable( TARGET  l95_addincrement.x
                        SOURCES AddIncrement.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_addincrement.x )

ecbuild_add_executable( TARGET  l95_eda.x
                        SOURCES EDA.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_eda.x )

ecbuild_add_executable( TARGET  l95_letkf.x
                        SOURCES LETKF.cc
                        LIBS    lorenz95)
oops_output_json_schema( l95_letkf.x )

ecbuild_add_executable( TARGET  l95_hofx.x
                        SOURCES HofX.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_hofx.x )

ecbuild_add_executable( TARGET  l95_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_hofx3d.x )

ecbuild_add_executable( TARGET  l95_enshofx.x
                        SOURCES EnsHofX.cc
                        LIBS    lorenz95 )
oops_output_json_schema( l95_enshofx.x )
