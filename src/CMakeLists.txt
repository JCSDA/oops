include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_src_files

oops/assimilation/BMatrix.h
oops/assimilation/CalcHofX.h
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
oops/assimilation/DRMinimizer.h
oops/assimilation/DRPCGMinimizer.h
oops/assimilation/DRPFOMMinimizer.h
oops/assimilation/DRPLanczosMinimizer.h
oops/assimilation/DualMinimizer.h
oops/assimilation/DualVector.h
oops/assimilation/FGMRES.h
oops/assimilation/FGMRESMinimizer.h
oops/assimilation/FtnTriDiagSpectrum.F90
oops/assimilation/FullGMRES.h
oops/assimilation/GMRESR.h
oops/assimilation/GMRESRMinimizer.h
oops/assimilation/HBHtMatrix.h
oops/assimilation/HessianMatrix.h
oops/assimilation/HMatrix.h
oops/assimilation/HtMatrix.h
oops/assimilation/HtRinvHMatrix.h
oops/assimilation/Increment4D.h
oops/assimilation/IncrementalAssimilation.h
oops/assimilation/instantiateCostFactory.h
oops/assimilation/instantiateMinFactory.h
oops/assimilation/IPCG.h
oops/assimilation/IPCGMinimizer.h
oops/assimilation/JqTerm.h
oops/assimilation/JqTermTLAD.h
oops/assimilation/LBGMRESRMinimizer.h
oops/assimilation/LBHessianMatrix.h
oops/assimilation/LBMinimizer.h
oops/assimilation/linsysteigen.h
oops/assimilation/Minimizer.h
oops/assimilation/MINRES.h
oops/assimilation/MINRESMinimizer.h
oops/assimilation/PCG.h
oops/assimilation/PCGMinimizer.h
oops/assimilation/PLanczos.h
oops/assimilation/PLanczosMinimizer.h
oops/assimilation/PrimalMinimizer.h
oops/assimilation/QNewtonLMP.h
oops/assimilation/RinvMatrix.h
oops/assimilation/rotmat.h
oops/assimilation/RPCGMinimizer.h
oops/assimilation/RPLanczosMinimizer.h
oops/assimilation/SaddlePointLMPMatrix.h
oops/assimilation/SaddlePointMatrix.h
oops/assimilation/SaddlePointMinimizer.h
oops/assimilation/SaddlePointPrecondMatrix.h
oops/assimilation/SaddlePointVector.h
oops/assimilation/SpectralLMP.h
oops/assimilation/State4D.h
oops/assimilation/TriDiagSolve.h
oops/assimilation/TriDiagSpectrum.h
oops/assimilation/UpHessSolve.h
oops/assimilation/UpTriSolve.h

oops/base/Accumulator.h
oops/base/Departures.h
oops/base/DeparturesEnsemble.h
oops/base/DiagonalMatrix.h
oops/base/DolphChebyshev.cc
oops/base/DolphChebyshev.h
oops/base/EnsembleCovariance.h
oops/base/EnsembleCovariance4D.h
oops/base/ErrorCovariance4D.h
oops/base/GeneralizedDepartures.h
oops/base/GeoVaLsWriter.h
oops/base/GridPoint.cc
oops/base/GridPoint.h
oops/base/HybridCovariance.h
oops/base/IdentityMatrix.h
oops/base/IncrementEnsemble.h
oops/base/instantiateCovarFactory.h
oops/base/instantiateObsFilterFactory.h
oops/base/LinearModelBase.h
oops/base/LinearVariableChangeBase.h
oops/base/Localization.h
oops/base/ModelBase.h
oops/base/ModelSpaceCovariance4DBase.h
oops/base/ModelSpaceCovarianceBase.h
oops/base/ObsAuxControls.h
oops/base/ObsAuxCovariances.h
oops/base/ObsAuxIncrements.h
oops/base/ObsEnsemble.h
oops/base/ObsErrorBase.h
oops/base/ObsErrors.h
oops/base/Observations.h
oops/base/Observer.h
oops/base/Observers.h
oops/base/ObserversTLAD.h
oops/base/ObserverTLAD.h
oops/base/ObsFilterBase.h
oops/base/ObsFilters.h
oops/base/ObsLocalizationBase.h
oops/base/ObsSpaceBase.cc
oops/base/ObsSpaceBase.h
oops/base/ObsSpaces.h
oops/base/PostBase.h
oops/base/PostBaseTLAD.h
oops/base/PostProcessor.h
oops/base/PostProcessorTLAD.h
oops/base/PostTimer.cc
oops/base/PostTimer.h
oops/base/PostTimerParameters.h
oops/base/StateEnsemble.h
oops/base/StateInfo.h
oops/base/StateWriter.h
oops/base/TrajectorySaver.h
oops/base/VariableChangeBase.h
oops/base/variables_f.cc
oops/base/variables_f.h
oops/base/variables_mod.F90
oops/base/Variables.cc
oops/base/Variables.h
oops/base/WeightedDiff.h
oops/base/WeightedDiffTLAD.h
oops/base/WeightedMean.h
oops/base/WeightingFct.cc
oops/base/WeightingFct.h

oops/generic/dcmip_initial_conditions_test_1_2_3_v5.f90
oops/generic/dcmip_initial_conditions_test_4_v3.f90
oops/generic/gc99.h
oops/generic/gc99.cc
oops/generic/IdLinearVariableChange.h
oops/generic/IdVariableChange.h
oops/generic/instantiateLocalizationFactory.h
oops/generic/instantiateModelFactory.h
oops/generic/instantiateObsErrorFactory.h
oops/generic/instantiateTlmFactory.h
oops/generic/instantiateVariableChangeFactory.h
oops/generic/LinearModelId.h
oops/generic/LocalizationGeneric.h
oops/generic/LocalObsErrorDiag.h
oops/generic/ObsErrorDiag.h
oops/generic/PseudoModel.h
oops/generic/unstructured_grid_f.h
oops/generic/unstructured_grid_interface.F90
oops/generic/unstructured_grid_mod.F90
oops/generic/UnstructuredGrid.cc
oops/generic/UnstructuredGrid.h

oops/interface/ErrorCovariance.h
oops/interface/Geometry.h
oops/interface/GeometryIterator.h
oops/interface/GeoVaLs.h
oops/interface/GetValues.h
oops/interface/Increment.h
oops/interface/LinearGetValues.h
oops/interface/LinearModel.h
oops/interface/LinearObsOperator.h
oops/interface/LinearVariableChange.h
oops/interface/LocalizationBase.h
oops/interface/Locations.h
oops/interface/Model.h
oops/interface/ModelAuxControl.h
oops/interface/ModelAuxCovariance.h
oops/interface/ModelAuxIncrement.h
oops/interface/ObsAuxControl.h
oops/interface/ObsAuxCovariance.h
oops/interface/ObsAuxIncrement.h
oops/interface/ObsDataVector.h
oops/interface/ObsDiagnostics.h
oops/interface/ObsErrorCovariance.h
oops/interface/ObsFilter.h
oops/interface/ObsLocalization.h
oops/interface/ObsOperator.h
oops/interface/ObsSpace.h
oops/interface/ObsVector.h
oops/interface/State.h
oops/interface/VariableChange.h

oops/parallel/mpi/mpi.cc
oops/parallel/mpi/mpi.h

oops/runs/AddIncrement.h
oops/runs/Application.h
oops/runs/ConvertState.h
oops/runs/DiffStates.h
oops/runs/Dirac.h
oops/runs/EDA.h
oops/runs/EnsForecasts.h
oops/runs/EnsHofX.h
oops/runs/EnsRecenter.h
oops/runs/EnsVariance.h
oops/runs/ExternalDFI.h
oops/runs/FindLocalObs.h
oops/runs/Forecast.h
oops/runs/GenEnsPertB.h
oops/runs/HofX.h
oops/runs/HofXNoModel.h
oops/runs/LETKF.h
oops/runs/LocalHofX.h
oops/runs/MakeObs.h
oops/runs/Run.cc
oops/runs/Run.h
oops/runs/StaticBInit.h
oops/runs/Test.h
oops/runs/Variational.h

oops/util/abor1_cpp.cc
oops/util/abor1_cpp.h
oops/util/abor1_ftn.F90
oops/util/config_f.cc
oops/util/config_f.h
oops/util/config_mod.F90
oops/util/config.intfb.h
oops/util/dateFunctions.cc
oops/util/dateFunctions.h
oops/util/datetime_f.cc
oops/util/datetime_f.h
oops/util/datetime_mod.F90
oops/util/DateTime.cc
oops/util/DateTime.h
oops/util/datetime.intfb.h
oops/util/dot_product.h
oops/util/duration_f.cc
oops/util/duration_f.h
oops/util/duration_mod.F90
oops/util/Duration.cc
oops/util/Duration.h
oops/util/duration.intfb.h
oops/util/Expect.h
oops/util/FloatCompare.h
oops/util/formats.h
oops/util/IntSetParser.cc
oops/util/IntSetParser.h
oops/util/IsAnyPointInVolumeInterior.h
oops/util/IsPointInVolumeInterior.h
oops/util/kinds.F90
oops/util/liboops_f.cc
oops/util/liboops_f.h
oops/util/liboops_mod.F90
oops/util/LibOOPS.cc
oops/util/LibOOPS.h
oops/util/Logger.h
oops/util/missing_values_f.cc
oops/util/missing_values_f.h
oops/util/missing_values_mod.F90
oops/util/netcdf_utils_mod.f90
oops/util/missingValues.cc
oops/util/missingValues.h
oops/util/ObjectCounter.h
oops/util/ObjectCountHelper.cc
oops/util/ObjectCountHelper.h
oops/util/Printable.h
oops/util/PrintAdjTest.h
oops/util/random_f.cc
oops/util/random_f.h
oops/util/random_mod.F90
oops/util/Random.h
oops/util/random.intfb.h
oops/util/ScalarOrMap.h
oops/util/signal_trap.cc
oops/util/sqr.h
oops/util/string_f_c_interface.cc
oops/util/string_f_c_interface.h
oops/util/string_f_c_mod.F90
oops/util/string_utils.F90
oops/util/stringFunctions.cc
oops/util/stringFunctions.h
oops/util/Timer.cc
oops/util/Timer.h
oops/util/TimerHelper.cc
oops/util/TimerHelper.h
oops/util/utilNamespaceDoc.h
oops/util/unstructured_interpolation_mod.f90

oops/util/parameters/OptionalParameter.h
oops/util/parameters/Parameter.h
oops/util/parameters/ParameterBase.cc
oops/util/parameters/ParameterBase.h
oops/util/parameters/Parameters.cc
oops/util/parameters/Parameters.h
oops/util/parameters/ParameterTraits.h
oops/util/parameters/ParameterTraitsScalarOrMap.h
oops/util/parameters/RequiredParameter.h

test/TestEnvironment.h
test/TestFixture.h

test/interface/ErrorCovariance.h
test/interface/Geometry.h
test/interface/GeometryIterator.h
test/interface/GeoVaLs.h
test/interface/GetValues.h
test/interface/Increment.h
test/interface/LinearGetValues.h
test/interface/LinearModel.h
test/interface/LinearObsOperator.h
test/interface/LinearVariableChange.h
test/interface/LocalObsSpace.h
test/interface/Locations.h
test/interface/Model.h
test/interface/ModelAuxControl.h
test/interface/ModelAuxCovariance.h
test/interface/ModelAuxIncrement.h
test/interface/ObsAuxControl.h
test/interface/ObsAuxCovariance.h
test/interface/ObsAuxIncrement.h
test/interface/ObsOperator.h
test/interface/ObsSpace.h
test/interface/ObsTestsFixture.h
test/interface/ObsVector.h
test/interface/State.h
test/interface/VariableChange.h

test/util/Fortran.h
test/util/IsAnyPointInVolumeInterior.h
test/util/Random.h
test/util/random.F90
test/util/FCString.h
test/util/f_c_string.F90
test/util/Parameters.h
test/util/ScalarOrMap.h
test/util/FloatCompare.h
)

list (APPEND oops_fheader_files
oops/util/linkedList_c.f
oops/util/linkedList_i.f
)

# macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${FILENAME} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)

# Link the oops test input files
list( APPEND oops_test_input
  test/testinput/random.yaml
  test/testinput/pushstringvector.yaml
  test/testinput/interp_interface.yaml
  test/testinput/variables.yaml
  test/testinput/parameters.yaml
  test/testinput/empty.yaml
)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/testinput)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${oops_test_input} )

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
                     LIBS       ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} Eigen3::Eigen fckit eckit_mpi eckit atlas atlas_f
                     INSTALL_HEADERS LISTED
                     TEMPLATES ${oops_fheader_files}
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                   )

ecbuild_add_test( TARGET oops_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

ecbuild_add_test( TARGET  test_oops_base_variables
                  SOURCES test/base/Variables.cc
                  ARGS    "test/testinput/variables.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_base_posttimer
                  SOURCES test/base/PostTimer.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops eckit)

ecbuild_add_test( TARGET  test_util_random
                  SOURCES test/base/Random.cc
                  ARGS    "test/testinput/random.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_pushstringvector
                  SOURCES test/base/FCString.cc
                  ARGS    "test/testinput/pushstringvector.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_parameters
                  SOURCES test/base/Parameters.cc
                  ARGS    "test/testinput/parameters.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_generic_gc99
                  SOURCES test/generic/gc99.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_isanypointinvolumeinterior
                  SOURCES test/util/IsAnyPointInVolumeInterior.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_datetime
                  SOURCES test/util/DateTime.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_duration
                  SOURCES test/util/Duration.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_intset_parser
                  SOURCES test/util/IntSetParser.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_scalarormap
                  SOURCES test/util/ScalarOrMap.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_floatcompare
                  SOURCES test/util/FloatCompare.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

if( HAVE_FCTEST )

add_fctest( TARGET  test_util_datetime_intfb.x
            SOURCES test/util/datetime.F90
            LINKER_LANGUAGE Fortran
            CONDITION HAVE_FCTEST
            LIBS    oops )

add_fctest( TARGET  test_util_duration_intfb.x
            SOURCES test/util/duration.F90
            LINKER_LANGUAGE Fortran
            CONDITION HAVE_FCTEST
            LIBS    oops )

endif()
