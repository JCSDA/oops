include(oops_output_json_schema)

ecbuild_add_executable( TARGET  qg_dfi.x
                        SOURCES qgDFI.cc
                        LIBS    qg
                       )
oops_output_json_schema( qg_dfi.x )

ecbuild_add_executable( TARGET  qg_ens_hofx.x
                        SOURCES qgEnsHofX.cc
                        LIBS    qg
                       )
oops_output_json_schema( qg_ens_hofx.x )

ecbuild_add_executable( TARGET  qg_ens_forecast.x
                        SOURCES qgEnsForecast.cc
                        LIBS    qg
                       )
oops_output_json_schema( qg_ens_forecast.x )

ecbuild_add_executable( TARGET  qg_forecast.x
                        SOURCES qgForecast.cc
                        LIBS    qg
                       )
oops_output_json_schema( qg_forecast.x )

ecbuild_add_executable( TARGET  qg_gen_ens_pert_B.x
                        SOURCES qgGenEnsPertB.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_gen_ens_pert_B.x )

ecbuild_add_executable( TARGET  qg_4dvar.x
                        SOURCES qg4DVar.cc
                        LIBS    qg
                      )
# TODO: oops_output_json_schema( qg_4dvar.x )

ecbuild_add_executable( TARGET  qg_hofx.x
                        SOURCES qgHofX.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_hofx.x )

ecbuild_add_executable( TARGET  qg_hofx3d.x
                        SOURCES qgHofX3D.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_hofx3d.x )

ecbuild_add_executable( TARGET  qg_diffstates.x
                        SOURCES qgDiffStates.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_diffstates.x )

ecbuild_add_executable( TARGET  qg_addincrement.x
                        SOURCES qgAddIncrement.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_addincrement.x )

ecbuild_add_executable( TARGET  qg_dirac.x
                        SOURCES qgDirac.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_dirac.x )

ecbuild_add_executable( TARGET  qg_staticbinit.x
                        SOURCES qgStaticBInit.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_staticbinit.x )

ecbuild_add_executable( TARGET  qg_eda.x
                        SOURCES qgEDA.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_eda.x )

ecbuild_add_executable( TARGET  qg_convertstate.x
                        SOURCES qgConvertState.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_convertstate.x )

ecbuild_add_executable( TARGET  qg_convertincrement.x
                        SOURCES qgConvertIncrement.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_convertincrement.x )

ecbuild_add_executable( TARGET  qg_ens_variance.x
                        SOURCES qgEnsVariance.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_ens_variance.x )

ecbuild_add_executable( TARGET  qg_letkf.x
                        SOURCES qgLETKF.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_letkf.x )

ecbuild_add_executable( TARGET  qg_ens_recenter.x
                        SOURCES qgEnsRecenter.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_ens_recenter.x )

ecbuild_add_executable( TARGET  qg_hybridgain.x
                        SOURCES qgHybridGain.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_hybridgain.x )

ecbuild_add_executable( TARGET  qg_rtpp.x
                        SOURCES qgRTPP.cc
                        LIBS    qg
                      )
oops_output_json_schema( qg_rtpp.x )
