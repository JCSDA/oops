ecbuild_add_executable( TARGET  qg_dfi.x
                        SOURCES qgDFI.cc
                        LIBS    qg
                       )

ecbuild_add_executable( TARGET  qg_ens_hofx.x
                        SOURCES qgEnsHofX.cc
                        LIBS    qg
                       )

ecbuild_add_executable( TARGET  qg_ens_forecast.x
                        SOURCES qgEnsForecast.cc
                        LIBS    qg
                       )

ecbuild_add_executable( TARGET  qg_forecast.x
                        SOURCES qgForecast.cc
                        LIBS    qg
                       )

ecbuild_add_executable( TARGET  qg_gen_ens_pert_B.x
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

ecbuild_add_executable( TARGET  qg_hofx_nomodel.x
                        SOURCES qgHofXNoModel.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_diffstates.x
                        SOURCES qgDiffStates.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_dirac.x
                        SOURCES qgDirac.cc
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

ecbuild_add_executable( TARGET  qg_ens_variance.x
                        SOURCES qgEnsVariance.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_letkf.x
                        SOURCES qgLETKF.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_ens_recenter.x
                        SOURCES qgEnsRecenter.cc
                        LIBS    qg
                      )

ecbuild_add_executable( TARGET  qg_rtpp.x
                        SOURCES qgRTPP.cc
                        LIBS    qg
                      )
