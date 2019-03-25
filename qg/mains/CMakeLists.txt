ecbuild_add_executable( TARGET  qg_dfi.x
                        SOURCES qgDFI.cc
                        LIBS    qg
                       )

ecbuild_add_executable( TARGET  qg_enshofx.x
                        SOURCES qgEnsHofX.cc
                        LIBS    qg
                       )

ecbuild_add_executable( TARGET  qg_forecast.x
                        SOURCES qgForecast.cc
                        LIBS    qg
                       )

ecbuild_add_executable( TARGET  qg_genpert.x
                        SOURCES qgGenEnsPertB.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_4dvar.x
                        SOURCES qg4DVar.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_hofx.x
                        SOURCES qgHofX.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_hofx3d.x
                        SOURCES qgHofX3D.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_makeobs.x
                        SOURCES qgMakeObs.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_dirac.x
                        SOURCES qgDirac.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_parameters.x
                        SOURCES qgEstimateParams.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_staticbinit.x
                        SOURCES qgStaticBInit.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_eda.x
                        SOURCES qgEDA.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_convertstate.x
                        SOURCES qgConvertState.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_ensvariance.x
                        SOURCES qgEnsVariance.cc
                        LIBS    qg
			)

		      
