ecbuild_add_executable( TARGET  qg_dfi.x
                        SOURCES qgDFI.cc
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

ecbuild_add_executable( TARGET  qg_makeobs.x
                        SOURCES qgMakeObs.cc
                        LIBS    qg
                      )

