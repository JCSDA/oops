include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_src_files

util/abor1_cpp.h
util/abor1_cpp.cc
util/abor1_ftn.F90
util/config.intfb.h
util/config_f.h
util/config_f.cc
util/config_mod.F90
util/dateFunctions.h
util/dateFunctions.cc
util/DateTime.h
util/DateTime.cc
util/datetime_f.h
util/datetime_f.cc
util/datetime.intfb.h
util/datetime_mod.F90
util/dot_product.h
util/Duration.h
util/Duration.cc
util/duration.intfb.h
util/duration_f.h
util/duration_f.cc
util/duration_mod.F90
util/formats.h
util/kinds.F90
util/LibOOPS.h
util/LibOOPS.cc
util/ObjectCounter.h
util/ObjectCountHelper.h
util/ObjectCountHelper.cc
util/Printable.h
util/PrintAdjTest.h
util/set_err_trap.F90
util/signal_trap.c
util/string_f_c_mod.F90
util/Timer.h
util/Timer.cc
util/TimerHelper.h
util/TimerHelper.cc
util/utilNamespaceDoc.h
util/Logger.h

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
oops/assimilation/JqTermAD.h
oops/assimilation/JqTermTL.h
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
oops/base/GeneralizedDepartures.h
oops/base/HybridCovariance.h
oops/base/IdentityMatrix.h
oops/base/LinearObsOperators.h
oops/base/ModelSpaceCovarianceBase.h
oops/base/ObsErrors.h
oops/base/ObsOperators.h
oops/base/ObsSpaces.h
oops/base/Observations.h
oops/base/Observer.h
oops/base/ObserverTL.h
oops/base/ObserverAD.h
oops/base/PostBase.h
oops/base/PostBaseTL.h
oops/base/PostBaseAD.h
oops/base/PostProcessor.h
oops/base/PostProcessorTL.h
oops/base/PostProcessorAD.h
oops/base/PostTimer.h
oops/base/PostTimer.cc
oops/base/StateInfo.h
oops/base/StateWriter.h
oops/base/TrajectorySaver.h
oops/base/WeightedDiff.h
oops/base/WeightedDiffTL.h
oops/base/WeightedDiffAD.h
oops/base/WeightedMean.h
oops/base/WeightingFct.h
oops/base/WeightingFct.cc
oops/base/instantiateCovarFactory.h

oops/generic/ObsErrorDiag.h
oops/generic/instantiateObsErrorFactory.h

oops/interface/ErrorCovariance.h
oops/interface/Geometry.h
oops/interface/Increment.h
oops/interface/LinearModel.h
oops/interface/LinearObsOperBase.h
oops/interface/LinearObsOperator.h
oops/interface/Localization.h
oops/interface/Locations.h
oops/interface/Model.h
oops/interface/ModelAtLocations.h
oops/interface/ModelAuxControl.h
oops/interface/ModelAuxCovariance.h
oops/interface/ModelAuxIncrement.h
oops/interface/ObsAuxControl.h
oops/interface/ObsAuxCovariance.h
oops/interface/ObsAuxIncrement.h
oops/interface/ObsErrorBase.h
oops/interface/ObsErrorCovariance.h
oops/interface/ObsOperator.h
oops/interface/ObsOperatorBase.h
oops/interface/ObsVector.h
oops/interface/ObservationSpace.h
oops/interface/State.h
oops/interface/Variables.h

oops/runs/Application.h
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

test/base/TestSuiteOpObsFixture.h
test/base/TestSuiteOpObsTLAD.h
test/base/TestSuiteOpObsTL.h
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
                     LIBS       ${LIBXML2_LIBRARIES} ${LAPACK_LIBRARIES} fckit eckit_mpi eckit
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                   )

ecbuild_add_test( TARGET  test_util.x
                  BOOST
                  SOURCES util/test/DateTime.cc
                          util/test/Duration.cc
                          util/test/test_main.cc
                  LIBS    oops )

if( HAVE_FCTEST )
add_fctest( TARGET  test_util_config_intfb.x
            SOURCES util/test/config.F90
            CONDITION HAVE_FCTEST
            LIBS    oops )

add_fctest( TARGET  test_util_datetime_intfb.x
            SOURCES util/test/datetime.F90
            CONDITION HAVE_FCTEST
            LIBS    oops )

add_fctest( TARGET  test_util_duration_intfb.x
            SOURCES util/test/duration.F90
            CONDITION HAVE_FCTEST
            LIBS    oops )
endif()

