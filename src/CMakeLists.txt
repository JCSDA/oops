include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_src_files

oops/assimilation/BMatrix.h
oops/assimilation/ControlIncrement.h
oops/assimilation/ControlVariable.h
oops/assimilation/CostFct3DVar.h
oops/assimilation/CostFct4DEnsVar.h
oops/assimilation/CostFct4DVar.h
oops/assimilation/CostFctWeak.h
oops/assimilation/CostFunction.h
oops/assimilation/CostJb3D.h
oops/assimilation/CostJb4D.h
oops/assimilation/CostJbJq.h
oops/assimilation/CostJbState.h
oops/assimilation/CostJbTotal.h
oops/assimilation/CostJcDFI.h
oops/assimilation/CostJo.h
oops/assimilation/CostTermBase.h
oops/assimilation/DRGMRESRMinimizer.h
oops/assimilation/DRIPCGMinimizer.h
oops/assimilation/DRPCGMinimizer.h
oops/assimilation/DRPLanczosMinimizer.h
oops/assimilation/DRMinimizer.h
oops/assimilation/DualMinimizer.h
oops/assimilation/DualVector.h
oops/assimilation/FGMRES.h
oops/assimilation/FGMRESMinimizer.h
oops/assimilation/FtnTriDiagSpectrum.F90
oops/assimilation/FullGMRES.h
oops/assimilation/GMRESR.h
oops/assimilation/GMRESRMinimizer.h
oops/assimilation/HBHtMatrix.h
oops/assimilation/LBGMRESRMinimizer.h
oops/assimilation/LBHessianMatrix.h
oops/assimilation/LBMinimizer.h
oops/assimilation/BMatrix.h
oops/assimilation/HMatrix.h
oops/assimilation/HessianMatrix.h
oops/assimilation/HtMatrix.h
oops/assimilation/HtRinvHMatrix.h
oops/assimilation/IPCG.h
oops/assimilation/IPCGMinimizer.h
oops/assimilation/Increment4D.h
oops/assimilation/IncrementalAssimilation.h
oops/assimilation/JqTerm.h
oops/assimilation/JqTermTLAD.h
oops/assimilation/MINRES.h
oops/assimilation/MINRESMinimizer.h
oops/assimilation/Minimizer.h
oops/assimilation/PCG.h
oops/assimilation/PCGMinimizer.h
oops/assimilation/PLanczos.h
oops/assimilation/PLanczosMinimizer.h
oops/assimilation/PrimalMinimizer.h
oops/assimilation/QNewtonLMP.h
oops/assimilation/RPCGMinimizer.h
oops/assimilation/RPLanczosMinimizer.h
oops/assimilation/RinvMatrix.h
oops/assimilation/SaddlePointLMPMatrix.h
oops/assimilation/SaddlePointMatrix.h
oops/assimilation/SaddlePointMinimizer.h
oops/assimilation/SaddlePointPrecondMatrix.h
oops/assimilation/SaddlePointVector.h
oops/assimilation/SpectralLMP.h
oops/assimilation/State4D.h
oops/assimilation/TriDiagSolve.h
oops/assimilation/TriDiagSpectrum.h
oops/assimilation/UpTriSolve.h
oops/assimilation/instantiateCostFactory.h
oops/assimilation/instantiateMinFactory.h
oops/assimilation/rotmat.h

oops/base/Accumulator.h
oops/base/Departures.h
oops/base/DiagonalMatrix.h
oops/base/DolphChebyshev.h
oops/base/DolphChebyshev.cc
oops/base/Ensemble.h
oops/base/EnsembleCovariance.h
oops/base/EnsemblesCollection.h
oops/base/FilterBase.h
oops/base/GeneralizedDepartures.h
oops/base/GeoVaLsWriter.h
oops/base/HybridCovariance.h
oops/base/IdentityMatrix.h
oops/base/InterpolatorTLAD.h
oops/base/LinearObsOperators.h
oops/base/ModelSpaceCovarianceBase.h
oops/base/ObsErrors.h
oops/base/ObsFilter.h
oops/base/ObsFilters.h
oops/base/ObsOperators.h
oops/base/ObsSpaces.h
oops/base/Observations.h
oops/base/Observer.h
oops/base/ObserverTLAD.h
oops/base/PostBase.h
oops/base/PostBaseTLAD.h
oops/base/PostProcessor.h
oops/base/PostProcessorTLAD.h
oops/base/PostTimer.h
oops/base/PostTimer.cc
oops/base/StateInfo.h
oops/base/StateWriter.h
oops/base/TrajectorySaver.h
oops/base/VariableChangeBase.h
oops/base/Variables.h
oops/base/Variables.cc
oops/base/WeightedDiff.h
oops/base/WeightedDiffTLAD.h
oops/base/WeightedMean.h
oops/base/WeightingFct.h
oops/base/WeightingFct.cc
oops/base/instantiateCovarFactory.h
oops/base/instantiateFilterFactory.h

oops/generic/ErrorCovarianceBUMP.h
oops/generic/LocalizationBUMP.h
oops/generic/ObsErrorDiag.h
oops/generic/UnstructuredGrid.cc
oops/generic/UnstructuredGrid.h
oops/generic/bump_f.h
oops/generic/bump_mod.F90
oops/generic/dcmip_initial_conditions_test_1_2_3_v5.f90
oops/generic/dcmip_initial_conditions_test_4_v3.f90
oops/generic/instantiateLocalizationFactory.h
oops/generic/instantiateObsErrorFactory.h
oops/generic/instantiateVariableChangeFactory.h
oops/generic/unstructured_grid_f.h
oops/generic/unstructured_grid_mod.F90
oops/generic/StatsVariableChange.h

oops/generic/bump/tools_const.F90
oops/generic/bump/tools_fit.F90
oops/generic/bump/tools_func.F90
oops/generic/bump/tools_icos.F90
oops/generic/bump/tools_kinds.F90
oops/generic/bump/tools_missing.F90
oops/generic/bump/tools_nc.F90
oops/generic/bump/tools_test.F90
oops/generic/bump/type_avg_blk.F90
oops/generic/bump/type_avg.F90
oops/generic/bump/type_bpar.F90
oops/generic/bump/type_bump.F90
oops/generic/bump/type_cmat_blk.F90
oops/generic/bump/type_cmat.F90
oops/generic/bump/type_com.F90
oops/generic/bump/type_cv_blk.F90
oops/generic/bump/type_cv.F90
oops/generic/bump/type_diag_blk.F90
oops/generic/bump/type_diag.F90
oops/generic/bump/type_displ.F90
oops/generic/bump/type_ens.F90
oops/generic/bump/type_geom.F90
oops/generic/bump/type_hdata.F90
oops/generic/bump/type_io.F90
oops/generic/bump/type_kdtree.F90
oops/generic/bump/type_lct_blk.F90
oops/generic/bump/type_lct.F90
oops/generic/bump/type_linop.F90
oops/generic/bump/type_mesh.F90
oops/generic/bump/type_minim.F90
oops/generic/bump/type_mom_blk.F90
oops/generic/bump/type_mom.F90
oops/generic/bump/type_mpl.F90
oops/generic/bump/type_nam.F90
oops/generic/bump/type_nicas_blk.F90
oops/generic/bump/type_nicas.F90
oops/generic/bump/type_obsop.F90
oops/generic/bump/type_rng.F90
oops/generic/bump/type_timer.F90
oops/generic/bump/external/tools_kdtree2.F90
oops/generic/bump/external/tools_kdtree2_pq.F90
oops/generic/bump/external/tools_qsort.F90
oops/generic/bump/external/tools_stripack.F90

oops/interface/ErrorCovariance.h
oops/interface/Geometry.h
oops/interface/GeoVaLs.h
oops/interface/Increment.h
oops/interface/InterpolatorTraj.h
oops/interface/LinearModel.h
oops/interface/LinearObsOperator.h
oops/interface/Localization.h
oops/interface/LocalizationBase.h
oops/interface/Locations.h
oops/interface/Model.h
oops/interface/ModelAuxControl.h
oops/interface/ModelAuxCovariance.h
oops/interface/ModelAuxIncrement.h
oops/interface/ObsAuxControl.h
oops/interface/ObsAuxCovariance.h
oops/interface/ObsAuxIncrement.h
oops/interface/ObsErrorBase.h
oops/interface/ObsErrorCovariance.h
oops/interface/ObsOperator.h
oops/interface/ObsSpaceBase.h
oops/interface/ObsVector.h
oops/interface/ObservationSpace.h
oops/interface/State.h
oops/interface/VariableChange.h

oops/runs/Application.h
oops/runs/Dirac.h
oops/runs/EnsForecasts.h
oops/runs/ExternalDFI.h
oops/runs/Forecast.h
oops/runs/GenEnsPertB.h
oops/runs/HofX.h
oops/runs/MakeObs.h
oops/runs/Run.h
oops/runs/Run.cc
oops/runs/Test.h
oops/runs/Variational.h

oops/util/abor1_cpp.h
oops/util/abor1_cpp.cc
oops/util/abor1_ftn.F90
oops/util/config.intfb.h
oops/util/config_f.h
oops/util/config_f.cc
oops/util/config_mod.F90
oops/util/dateFunctions.h
oops/util/dateFunctions.cc
oops/util/DateTime.h
oops/util/DateTime.cc
oops/util/datetime_f.h
oops/util/datetime_f.cc
oops/util/datetime.intfb.h
oops/util/datetime_mod.F90
oops/util/dot_product.h
oops/util/Duration.h
oops/util/Duration.cc
oops/util/duration.intfb.h
oops/util/duration_f.h
oops/util/duration_f.cc
oops/util/duration_mod.F90
oops/util/formats.h
oops/util/kinds.F90
oops/util/LibOOPS.h
oops/util/LibOOPS.cc
oops/util/Logger.h
oops/util/ObjectCounter.h
oops/util/ObjectCountHelper.h
oops/util/ObjectCountHelper.cc
oops/util/Printable.h
oops/util/PrintAdjTest.h
oops/util/printStackTrace.h
oops/util/set_err_trap.F90
oops/util/signal_trap.c
oops/util/string_f_c_mod.F90
oops/util/Timer.h
oops/util/Timer.cc
oops/util/TimerHelper.h
oops/util/TimerHelper.cc
oops/util/utilNamespaceDoc.h

test/interface/GeoVaLs.h
test/interface/LinearObsOperator.h
test/interface/Locations.h
test/interface/ObservationSpace.h
test/interface/ObsOperator.h
test/interface/ObsTestsFixture.h
test/interface/ObsVector.h
test/interface/VariableChange.h
test/TestEnvironment.h

)

function(generate_fortran_bindings output filename)

  set( options "" )
  set( single_value_args OUTPUT MODULE )
  set( multi_value_args "" )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  get_filename_component(base ${filename} NAME_WE)
  set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
  set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}_c_binding.f90)

  if( _PAR_OUTPUT )
    set(outfile ${_PAR_OUTPUT})
  endif()
  set(${output} ${${output}} ${outfile} PARENT_SCOPE)

  if( _PAR_MODULE )
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND python ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile} -m ${_PAR_MODULE}
        DEPENDS ${filename} )
  else()
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND python ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
        DEPENDS ${filename} )
  endif()
  set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
endfunction()

ecbuild_add_library( TARGET     oops
                     SOURCES    
                                ${FORTRAN_BINDINGS}
                                ${oops_src_files}
                     LIBS       ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} fckit eckit_mpi eckit
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                   )

ecbuild_add_test( TARGET oops_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

