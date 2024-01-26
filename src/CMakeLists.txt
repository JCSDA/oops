include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_src_files

oops/assimilation/BMatrix.h
oops/assimilation/CMatrix.h
oops/assimilation/ControlIncrement.h
oops/assimilation/ControlVariable.h
oops/assimilation/CostFct3DVar.h
oops/assimilation/CostFct4DEnsVar.h
oops/assimilation/CostFct4DVar.h
oops/assimilation/CostFctFGAT.h
oops/assimilation/CostFctWeak.h
oops/assimilation/CostFunction.h
oops/assimilation/CostJb3D.h
oops/assimilation/CostJb4D.h
oops/assimilation/CostJbJq.h
oops/assimilation/CostJbModelAux.h
oops/assimilation/CostJbObsAux.h
oops/assimilation/CostJbState.h
oops/assimilation/CostJbTotal.h
oops/assimilation/CostJcDFI.h
oops/assimilation/CostJo.h
oops/assimilation/CostPert.h
oops/assimilation/CostTermBase.h
oops/assimilation/DRGMRESRMinimizer.h
oops/assimilation/DRIPCGMinimizer.h
oops/assimilation/DRMinimizer.h
oops/assimilation/DRPBlockLanczosMinimizer.h
oops/assimilation/DRPCGMinimizer.h
oops/assimilation/DRPFOMMinimizer.h
oops/assimilation/DRPLanczosMinimizer.h
oops/assimilation/DualMinimizer.h
oops/assimilation/DualVector.h
oops/assimilation/FGMRES.h
oops/assimilation/FGMRESMinimizer.h
oops/assimilation/FtnTriDiagSpectrum.F90
oops/assimilation/FullGMRES.h
oops/assimilation/GETKFSolver.h
oops/assimilation/gletkfInterface.h
oops/assimilation/gletkf_interface.f90
oops/assimilation/gletkf_mod.f90
oops/assimilation/GMRESR.h
oops/assimilation/GMRESRMinimizer.h
oops/assimilation/HBHtMatrix.h
oops/assimilation/HessianMatrix.h
oops/assimilation/HMatrix.h
oops/assimilation/HtMatrix.h
oops/assimilation/HtRinvHMatrix.h
oops/assimilation/IncrementalAssimilation.h
oops/assimilation/instantiateCostFactory.h
oops/assimilation/instantiateLocalEnsembleSolverFactory.h
oops/assimilation/instantiateMinFactory.h
oops/assimilation/IPCG.h
oops/assimilation/IPCGMinimizer.h
oops/assimilation/JqTerm.h
oops/assimilation/JqTermTLAD.h
oops/assimilation/LBGMRESRMinimizer.h
oops/assimilation/LBHessianMatrix.h
oops/assimilation/LBMinimizer.h
oops/assimilation/LETKFSolver.h
oops/assimilation/LETKFSolverGSI.h
oops/assimilation/LocalEnsembleSolverParameters.h
oops/assimilation/LocalEnsembleSolver.h
oops/assimilation/Minimizer.h
oops/assimilation/MinimizerUtils.cc
oops/assimilation/MinimizerUtils.h
oops/assimilation/MINRES.h
oops/assimilation/MINRESMinimizer.h
oops/assimilation/PCG.h
oops/assimilation/PCGMinimizer.h
oops/assimilation/PFF.h
oops/assimilation/PLanczos.h
oops/assimilation/PLanczosMinimizer.h
oops/assimilation/PMatrix.h
oops/assimilation/PrimalMinimizer.h
oops/assimilation/QNewtonLMP.h
oops/assimilation/RinvHMatrix.h
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
oops/assimilation/TriDiagSolve.h
oops/assimilation/TriDiagSpectrum.h
oops/assimilation/UpHessSolve.h
oops/assimilation/UpTriSolve.h
oops/base/Accumulator.h
oops/base/AnalyticInit.h
oops/base/DataSetBase.h
oops/base/Departures.h
oops/base/DeparturesEnsemble.h
oops/base/DiagonalMatrix.h
oops/base/DolphChebyshev.cc
oops/base/DolphChebyshev.h
oops/base/EnsembleCovariance.h
oops/base/FieldSet3D.cc
oops/base/FieldSet3D.h
oops/base/FieldSet4D.cc
oops/base/FieldSet4D.h
oops/base/FieldSets.cc
oops/base/FieldSets.h
oops/base/ForecastParameters.h
oops/base/GeneralizedDepartures.h
oops/base/Geometry.h
oops/base/GeometryData.cc
oops/base/GeometryData.h
oops/base/GetValuePosts.h
oops/base/GetValues.h
oops/base/GetValueTLADs.h
oops/base/LocalIncrement.cc
oops/base/LocalIncrement.h
oops/base/Locations.h
oops/base/HybridCovariance.h
oops/base/IdentityMatrix.h
oops/base/Increment.h
oops/base/Increment4D.h
oops/base/IncrementEnsemble.h
oops/base/IncrementEnsemble4D.h
oops/base/IncrementSet.h
oops/base/Inflation.h
oops/base/InflationBase.h
oops/base/instantiateCovarFactory.h
oops/base/instantiateInflationFactory.h
oops/base/instantiateObsFilterFactory.h
oops/base/LatLonGridPostProcessor.h
oops/base/LatLonGridWriter.h
oops/base/LinearModel.h
oops/base/LinearVariableChangeParametersBase.h
oops/base/Localization.h
oops/base/Model.h
oops/base/ModelSpaceCovarianceBase.h
oops/base/ModelSpaceCovarianceParametersBase.h
oops/base/MultiplicativeInflation.h
oops/base/NormBase.h
oops/base/ObsAuxControls.h
oops/base/ObsAuxCovariances.h
oops/base/ObsAuxIncrements.h
oops/base/ObsAuxPreconditioners.h
oops/base/ObsEnsemble.h
oops/base/ObsError.h
oops/base/ObsErrors.h
oops/base/Observations.h
oops/base/Observer.h
oops/base/Observers.h
oops/base/ObserversTLAD.h
oops/base/ObserverTLAD.h
oops/base/ObserverUtils.h
oops/base/ObsFilter.h
oops/base/ObsFilters.h
oops/base/ObsLocalizationBase.h
oops/base/ObsLocalizationParametersBase.h
oops/base/ObsLocalizations.h
oops/base/ObsOperatorBase.h
oops/base/ObsOperatorPert.h
oops/base/ObsSpaceBase.cc
oops/base/ObsSpaceBase.h
oops/base/ObsSpaces.h
oops/base/ObsTypeParameters.h
oops/base/ObsVector.h
oops/base/ParameterTraitsVariables.cc
oops/base/ParameterTraitsVariables.h
oops/base/PostBase.h
oops/base/PostBaseTLAD.h
oops/base/PostProcessor.h
oops/base/PostProcessorTLAD.h
oops/base/PostTimer.cc
oops/base/PostTimer.h
oops/base/PostTimerParameters.h
oops/base/RTPP.h
oops/base/RTPS.h
oops/base/SamplingMethodSelector.h
oops/base/State.h
oops/base/State4D.h
oops/base/StateEnsemble.h
oops/base/StateEnsemble4D.h
oops/base/StateInfo.h
oops/base/StateSaver.h
oops/base/StateSet.h
oops/base/StateWriter.h
oops/base/TrajectorySaver.h
oops/base/TrivialSamplingMethodSelector.h
oops/base/WriteParametersBase.cc
oops/base/WriteParametersBase.h
oops/base/VariableChangeParametersBase.h
oops/base/variables_f.cc
oops/base/variables_f.h
oops/base/variables_mod.F90
oops/base/Variables.cc
oops/base/Variables.h
oops/base/WeightedDiff.h
oops/base/WeightedDiffTLAD.h
oops/base/WeightedMean.h
oops/base/WeightingFct.h

oops/contrib/dcmip_initial_conditions_test_1_2_3_v5.f90
oops/contrib/dcmip_initial_conditions_test_4_v3.f90
oops/contrib/dcmip_initial_conditions_test_4_v3_interface.f90
oops/contrib/dcmipInitialConditionsTest_4_v3_Fortran.h

oops/coupled/AuxCoupledModel.h
oops/coupled/GeometryCoupled.h
oops/coupled/GetValuesCoupled.h
oops/coupled/IncrementCoupled.h
oops/coupled/ModelCoupled.h
oops/coupled/StateCoupled.h
oops/coupled/TraitCoupled.h
oops/coupled/UtilsCoupled.cc
oops/coupled/UtilsCoupled.h
oops/coupled/VariableChangeCoupled.h
oops/coupled/instantiateCoupledFactory.h

oops/external/stripack/stripack.cc
oops/external/stripack/stripack.h
oops/external/stripack/stripack.f90
oops/external/stripack/stripack_interface.f90

oops/generic/AnalyticInitBase.h
oops/generic/AtlasInterpolator.cc
oops/generic/AtlasInterpolator.h
oops/generic/fft_gpoint2spectral_f.F90
oops/generic/fft_init_f.F90
oops/generic/fft_interface_f.h
oops/generic/fft_multiple_f.F90
oops/generic/fft_pass_proc_f.F90
oops/generic/fft_spectral2gpoint_f.F90
oops/generic/gc99.h
oops/generic/gc99.cc
oops/generic/gc99_mod.F90
oops/generic/GeoVaLsWriter.h
oops/generic/GlobalInterpolator.cc
oops/generic/GlobalInterpolator.h
oops/generic/GlobalTemplatedInterpolator.h
oops/generic/HtlmCalculator.h
oops/generic/HtlmEnsemble.h
oops/generic/HtlmRegularization.cc
oops/generic/HtlmRegularization.h
oops/generic/HybridLinearModel.h
oops/generic/HybridLinearModelCoeffs.h
oops/generic/IdentityLinearModel.h
oops/generic/IdentityModel.h
oops/generic/LocalInterpolatorBase.cc
oops/generic/LocalInterpolatorBase.h
oops/generic/instantiateLinearModelFactory.h
oops/generic/instantiateModelFactory.h
oops/generic/instantiateNormFactory.h
oops/generic/instantiateObsErrorFactory.h
oops/generic/LinearModelBase.h
oops/generic/LocalizationBase.h
oops/generic/ModelBase.h
oops/generic/NormScalar.h
oops/generic/ObsErrorDiag.h
oops/generic/ObsErrorDiag.cc
oops/generic/ObsErrorBase.h
oops/generic/ObsFilterBase.h
oops/generic/ObsFilterParametersBase.h
oops/generic/PseudoModel.h
oops/generic/PseudoModelState4D.h
oops/generic/SimpleLinearModel.h
oops/generic/SimpleLinearModelMultiresolution.h
oops/generic/SimpleLinearModelResidualForm.h
oops/generic/soar.h
oops/generic/soar.cc
oops/generic/UnstructuredInterpolator.h
oops/generic/UnstructuredInterpolator.cc
oops/generic/VerticalLocEV.h

oops/interface/AnalyticInitBase.h
oops/interface/ErrorCovariance.h
oops/interface/Geometry.h
oops/interface/GeometryIterator.h
oops/interface/GeoVaLs.h
oops/interface/Increment.h
oops/interface/LinearModelBase.h
oops/interface/LinearObsOperator.h
oops/interface/LinearVariableChange.h
oops/interface/LocalInterpolator.h
oops/interface/LocalizationBase.h
oops/interface/ModelAuxControl.h
oops/interface/ModelAuxCovariance.h
oops/interface/ModelAuxIncrement.h
oops/interface/ModelBase.h
oops/interface/ModelData.h
oops/interface/NormGradient.h
oops/interface/ObsAuxControl.h
oops/interface/ObsAuxCovariance.h
oops/interface/ObsAuxIncrement.h
oops/interface/ObsAuxPreconditioner.h
oops/interface/ObsDataVector.h
oops/interface/ObsDataVector_head.h
oops/interface/ObsDiagnostics.h
oops/interface/ObsErrorBase.h
oops/interface/ObsFilterBase.h
oops/interface/ObsOperator.h
oops/interface/ObsSpace.h
oops/interface/ObsVector.h
oops/interface/SampledLocations.h
oops/interface/State.h
oops/interface/VariableChange.h

oops/mpi/mpi.cc
oops/mpi/mpi.h

oops/runs/AddIncrement.h
oops/runs/AdjointForecast.h
oops/runs/Application.cc
oops/runs/Application.h
oops/runs/ControlPert.h
oops/runs/ConvertIncrement.h
oops/runs/ConvertState.h
oops/runs/ConvertToLatLon.h
oops/runs/DiffStates.h
oops/runs/EnsembleApplication.h
oops/runs/EnsembleInflation.h
oops/runs/EnsMeanAndVariance.h
oops/runs/EnsRecenter.h
oops/runs/ExternalDFI.h
oops/runs/Forecast.h
oops/runs/GenEnsPertB.h
oops/runs/GenHybridLinearModelCoeffs.h
oops/runs/HofX3D.h
oops/runs/HofX4D.h
oops/runs/HybridGain.h
oops/runs/LinearizationError.h
oops/runs/LocalEnsembleDA.h
oops/runs/Run.cc
oops/runs/Run.h
oops/runs/SqrtOfVertLoc.h
oops/runs/Test.h
oops/runs/Variational.h

oops/util/abor1_cpp.cc
oops/util/abor1_cpp.h
oops/util/abor1_ftn.F90
oops/util/algorithms.h
oops/util/AnyOf.h
oops/util/AssociativeContainers.h
oops/util/AtlasArrayUtil.cc
oops/util/AtlasArrayUtil.h
oops/util/CompareNVectors.h
oops/util/CompositePath.cc
oops/util/CompositePath.h
oops/util/ConfigFunctions.cc
oops/util/ConfigFunctions.h
oops/util/ConfigHelpers.cc
oops/util/ConfigHelpers.h
oops/util/dateFunctions.cc
oops/util/dateFunctions.h
oops/util/datetime_f.cc
oops/util/datetime_f.h
oops/util/datetime_mod.F90
oops/util/PartialDateTime.cc
oops/util/PartialDateTime.h
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
oops/util/FieldSetHelpers.cc
oops/util/FieldSetHelpers.h
oops/util/FieldSetOperations.cc
oops/util/FieldSetOperations.h
oops/util/FloatCompare.h
oops/util/FunctionSpaceHelpers.cc
oops/util/FunctionSpaceHelpers.h
oops/util/formats.h
oops/util/gatherPrint.h
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
oops/util/LocalEnvironment.cc
oops/util/LocalEnvironment.h
oops/util/Logger_f.cc
oops/util/Logger_f.h
oops/util/Logger.h
oops/util/logger_mod.F90
oops/util/missing_values_f.cc
oops/util/missing_values_f.h
oops/util/missing_values_mod.F90
oops/util/missingValues.cc
oops/util/missingValues.h
oops/util/NamedEnumerator.h
oops/util/netcdf_utils_mod.f90
oops/util/ObjectCounter.h
oops/util/ObjectCountHelper.cc
oops/util/ObjectCountHelper.h
oops/util/Printable.h
oops/util/PrintAdjTest.h
oops/util/printRunStats.cc
oops/util/printRunStats.h
oops/util/PropertiesOfNVectors.h
oops/util/random_f.cc
oops/util/random_f.h
oops/util/random_mod.F90
oops/util/Random.h
oops/util/random.intfb.h
oops/util/RandomField.h
oops/util/Range.h
oops/util/ScalarOrMap.h
oops/util/Serializable.h
oops/util/signal_trap.h
oops/util/signal_trap.cc
oops/util/sqr.h
oops/util/string_f_c_interface.cc
oops/util/string_f_c_interface.h
oops/util/string_f_c_mod.F90
oops/util/string_utils.F90
oops/util/stringFunctions.cc
oops/util/stringFunctions.h
oops/util/TestReference.cc
oops/util/TestReference.h
oops/util/Timer.cc
oops/util/Timer.h
oops/util/TimerHelper.cc
oops/util/TimerHelper.h
oops/util/TimeWindow.cc
oops/util/TimeWindow.h
oops/util/TypeTraits.h
oops/util/utilNamespaceDoc.h
oops/util/wildcard.h
oops/util/wildcard.cc
oops/util/workflow.h
oops/util/workflow.cc
oops/util/WorkflowUpdater.h

oops/util/parameters/ConfigurationParameter.cc
oops/util/parameters/ConfigurationParameter.h
oops/util/parameters/EmptyParameters.h
oops/util/parameters/GenericParameters.h
oops/util/parameters/HasDiracParameters_.h
oops/util/parameters/HasParameters_.h
oops/util/parameters/HasReadParameters_.h
oops/util/parameters/HasWriteParameters_.h
oops/util/parameters/IgnoreOtherParameters.cc
oops/util/parameters/IgnoreOtherParameters.h
oops/util/parameters/NumericConstraints.h
oops/util/parameters/ObjectJsonSchema.cc
oops/util/parameters/ObjectJsonSchema.h
oops/util/parameters/OptionalParameter.h
oops/util/parameters/OptionalParameter.cc
oops/util/parameters/OptionalPolymorphicParameter.h
oops/util/parameters/Parameter.h
oops/util/parameters/ParameterBase.cc
oops/util/parameters/ParameterBase.h
oops/util/parameters/ParameterConstraint.h
oops/util/parameters/Parameters.cc
oops/util/parameters/Parameters.h
oops/util/parameters/ParametersOrConfiguration.h
oops/util/parameters/ParameterTraits.cc
oops/util/parameters/ParameterTraits.h
oops/util/parameters/ParameterTraitsAnyOf.h
oops/util/parameters/ParameterTraitsScalarOrMap.h
oops/util/parameters/PolymorphicParameter.h
oops/util/parameters/PolymorphicParameterTraits.h
oops/util/parameters/PropertyJsonSchema.cc
oops/util/parameters/PropertyJsonSchema.h
oops/util/parameters/RequiredParameter.h
oops/util/parameters/RequiredPolymorphicParameter.h
)

list( APPEND oops_test_src_files
test/TestEnvironment.h
test/TestFixture.h

test/assimilation/FullGMRES.h
test/assimilation/rotmat.h
test/assimilation/SolveMatrixEquation.h
test/assimilation/SpectralLMP.h
test/assimilation/TestVector3D.h
test/assimilation/TriDiagSolve.h
test/assimilation/Vector3D.cc
test/assimilation/Vector3D.h

test/base/Fortran.h
test/base/NormBase.h
test/base/ObsErrorCovariance.h
test/base/ObsLocalizations.h
test/base/variables.F90

test/interface/ErrorCovariance.h
test/interface/FieldSet4D.h
test/interface/Geometry.h
test/interface/GeometryIterator.h
test/interface/GeometryFixture.h
test/interface/GeoVaLs.h
test/interface/GetValues.h
test/interface/GlobalInterpolator.h
test/interface/Increment.h
test/interface/LinearModel.h
test/interface/LinearObsOperator.h
test/interface/LinearVariableChange.h
test/interface/Model.h
test/interface/ModelAuxControl.h
test/interface/ModelAuxCovariance.h
test/interface/ModelAuxIncrement.h
test/interface/ObsAuxControl.h
test/interface/ObsAuxCovariance.h
test/interface/ObsAuxIncrement.h
test/interface/ObsDataVector.h
test/interface/ObsIterator.h
test/interface/ObsLocalization.h
test/interface/ObsOperator.h
test/interface/ObsSpace.h
test/interface/ObsTestsFixture.h
test/interface/ObsVector.h
test/interface/SampledLocations.h
test/interface/State.h
test/interface/UnstructuredInterpolator.h
test/interface/VariableChange.h

test/mpi/mpi.h

test/generic/fft_multiple.h
test/generic/PseudoModelState4D.h
test/generic/VerticalLocEV.h

test/util/ArrayUtil.h
test/util/Fortran.h
test/util/IsAnyPointInVolumeInterior.h
test/util/Random.h
test/util/random.F90
test/util/FCString.h
test/util/f_c_string.F90
test/util/AssociativeContainers.h
test/util/Parameters.h
test/util/ScalarOrMap.h
test/util/FloatCompare.h
test/util/CompareNVectors.h
test/util/CompositePath.h
test/util/MissingValues.h
test/util/PropertiesOfNVectors.h
test/util/stringFunctions.h
test/util/LocalEnvironment.h
test/util/TestReference.h
test/util/TypeTraits.h
test/util/algorithms.h
test/util/ConfigFunctions.h
test/util/ConfigHelpers.h
test/util/TimeWindow.h
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

# oops test input files
list( APPEND oops_test_input
  test/testinput/random.yaml
  test/testinput/pushstringvector.yaml
  test/testinput/variables.yaml
  test/testinput/parameters.yaml
  test/testinput/empty.yaml
  test/testinput/mpi.yaml
  test/testinput/spectrallmp.yaml
  test/testinput/fft_multiple.yaml
  test/testinput/hello.yaml
)
# oops test output files
list( APPEND oops_test_output
  test/testoutput/fft_multiple.test
)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/testinput)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${oops_test_input} )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/testoutput)
CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${oops_test_output} )

ecbuild_add_library( TARGET     ${PROJECT_NAME}
                     SOURCES    ${oops_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
                     LINKER_LANGUAGE CXX )

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/oops/util/linkedList_i.f
               ${CMAKE_CURRENT_SOURCE_DIR}/oops/util/linkedList_c.f
               DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oops/util )

target_sources( ${PROJECT_NAME} PRIVATE ${oops_test_src_files})
install( DIRECTORY test DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oops FILES_MATCHING PATTERN "*.h")

if(OpenMP_FOUND)
    target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX )
endif()

# Mark Boost as a system library to silence compiler warnings from Boost headers
target_include_directories( ${PROJECT_NAME} SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} )

target_link_libraries( ${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C ) # For utils/netcdf_utils_mod
target_link_libraries( ${PROJECT_NAME} PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} )
target_link_libraries( ${PROJECT_NAME} PUBLIC Eigen3::Eigen )
target_link_libraries( ${PROJECT_NAME} PUBLIC eckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC fckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC atlas_f )
target_link_libraries( ${PROJECT_NAME} PUBLIC Boost::boost )
if ( GPTL_FOUND )
	target_link_libraries( ${PROJECT_NAME} PUBLIC GPTL::GPTL )
endif()

# Stack traces on floating point exceptions
include( backtrace_deps )

target_link_libraries( ${PROJECT_NAME} PUBLIC ${OOPS_STACKTRACE_${OOPS_STACKTRACE_PROVIDER}_LIBS} )
target_compile_definitions( ${PROJECT_NAME} PRIVATE ${OOPS_STACKTRACE_${OOPS_STACKTRACE_PROVIDER}_DEFS} )

#Configure include directory layout for build-tree to match install-tree
set(OOPS_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(oops_test_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory ${OOPS_BUILD_DIR_INCLUDE_PATH}/oops
                                           COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/test ${OOPS_BUILD_DIR_INCLUDE_PATH}/oops/test)
add_dependencies(oops oops_test_headers)

target_include_directories(${PROJECT_NAME} INTERFACE
                                $<BUILD_INTERFACE:${OOPS_BUILD_DIR_INCLUDE_PATH}>
                                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                $<INSTALL_INTERFACE:include/oops>)

#Fortran module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                                $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                                $<INSTALL_INTERFACE:${MODULE_DIR}>)

if (${nlohmann_json_FOUND} AND ${nlohmann_json_schema_validator_FOUND})
  message(STATUS "[OOPS] JSON Schema Validator activated.")
  target_link_libraries( oops PRIVATE nlohmann_json nlohmann_json_schema_validator )
  target_compile_definitions( oops PRIVATE OOPS_HAVE_NLOHMANN_JSON_SCHEMA_VALIDATOR )
endif()

## Tests

ecbuild_add_test( TARGET oops_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

ecbuild_add_executable(
  TARGET  test_oops_base_dummy_run
  SOURCES test/base/DummyRun.cc
  LIBS oops
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_one
  COMMAND test_oops_base_dummy_run
  ARGS    "test/testinput/hello.yaml"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_one PROPERTIES
  PASS_REGULAR_EXPRESSION "hello earth"
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_no_validate
  COMMAND test_oops_base_dummy_run
  ARGS    "--no-validate" "test/testinput/hello.yaml"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_one PROPERTIES
  PASS_REGULAR_EXPRESSION "hello earth"
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_validate_zero
  COMMAND test_oops_base_dummy_run
  ARGS    "--validate-only" "test/testinput/hello.yaml"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_validate_zero PROPERTIES
  PASS_REGULAR_EXPRESSION "Configuration OK"
)

if (${nlohmann_json_FOUND} AND ${nlohmann_json_schema_validator_FOUND})
  ecbuild_add_test(
    TARGET  test_oops_base_dummy_run_validate_one
    COMMAND test_oops_base_dummy_run
    ARGS    "--validate-only" "test/testinput/empty.yaml"
    DEPENDS test_oops_base_dummy_run
  )
  set_tests_properties(test_oops_base_dummy_run_validate_one PROPERTIES
    WILL_FAIL true
  )
endif()

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_zero
  COMMAND test_oops_base_dummy_run
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_zero PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_one
  COMMAND test_oops_base_dummy_run
                  ARGS    "--validate-only"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_one PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_two
  COMMAND test_oops_base_dummy_run
  ARGS    "--output-json-schema=x" "test/testinput/hello.yaml"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_two PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_bad_arg_three
  COMMAND test_oops_base_dummy_run
  ARGS    "--validate-only" "test/testinput/hello.yaml" "out"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_bad_arg_three PROPERTIES
  WILL_FAIL true
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_help
  COMMAND test_oops_base_dummy_run
  ARGS    "--help"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_help PROPERTIES
  PASS_REGULAR_EXPRESSION "Usages:"
)

ecbuild_add_test(
  TARGET  test_oops_base_dummy_run_h
  COMMAND test_oops_base_dummy_run
  ARGS    "-h"
  DEPENDS test_oops_base_dummy_run
)
set_tests_properties(test_oops_base_dummy_run_h PROPERTIES
  PASS_REGULAR_EXPRESSION "Usages:"
)

ecbuild_add_test( TARGET  test_oops_base_variables
                  SOURCES test/base/Variables.cc
                  ARGS    "test/testinput/variables.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_base_posttimer
                  SOURCES test/base/PostTimer.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops eckit)

ecbuild_add_test( TARGET  test_util_signal_trap_fpe_div_by_zero
                  SOURCES test/util/signal_trap_fpe_div_by_zero.cc
                  LIBS    oops)
set_property( TEST test_util_signal_trap_fpe_div_by_zero PROPERTY WILL_FAIL TRUE )

ecbuild_add_test( TARGET  test_util_signal_trap_fpe_invalid_op
                  SOURCES test/util/signal_trap_fpe_invalid_op.cc
                  LIBS    oops)
set_property( TEST test_util_signal_trap_fpe_invalid_op PROPERTY WILL_FAIL TRUE )

# This test enables fpe trapping but doesn't perform any "bad" operations. So
# the expected outcome is for the test program to successfully complete.
ecbuild_add_test( TARGET  test_util_signal_trap_fpe_valid_op
                  SOURCES test/util/signal_trap_fpe_valid_op.cc
                  LIBS    oops)

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

ecbuild_add_test( TARGET  test_generic_soar
                  SOURCES test/generic/soar.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_coupled_splitvariables
                  SOURCES test/coupled/SplitVariables.cc
                  LIBS    oops eckit )

ecbuild_add_test( TARGET  test_util_isanypointinvolumeinterior
                  SOURCES test/util/IsAnyPointInVolumeInterior.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_partialdatetime
                  SOURCES test/util/PartialDateTime.cc
                  LIBS    oops eckit )

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

ecbuild_add_test( TARGET  test_util_compositepath
                  SOURCES test/util/CompositePath.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_stringfunctions
                  SOURCES test/util/stringFunctions.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_testreference
                  SOURCES test/util/TestReference.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_range
                  SOURCES test/util/Range.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_mpi_mpi
                  MPI     4
                  SOURCES test/mpi/mpi.cc
                  ARGS    "test/testinput/mpi.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_fft_multiple
                  SOURCES test/generic/fft_multiple.cc
                  ARGS    "test/testinput/fft_multiple.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_algorithms
                  SOURCES test/util/algorithms.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_comparenvectors
                  SOURCES test/util/CompareNVectors.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_missingvalues
                  SOURCES test/util/MissingValues.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_associativecontainers
                  SOURCES test/util/AssociativeContainers.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_propertiesofnvectors
                  SOURCES test/util/PropertiesOfNVectors.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_localenvironment
                  SOURCES test/util/LocalEnvironment.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_typetraits
                  SOURCES test/util/TypeTraits.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_wildcard
                  SOURCES test/util/wildcard.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_configfunctions
                  SOURCES test/util/ConfigFunctions.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_confighelpers
                  SOURCES test/util/ConfigHelpers.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_timewindow
                  SOURCES test/util/TimeWindow.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_arrayutil
                  MPI     2
                  SOURCES test/util/ArrayUtil.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_fieldset_helpers_and_operations
                  SOURCES test/util/FieldSetHelpersAndOperations.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_functionspace_helpers
                  SOURCES test/util/FunctionSpaceHelpers.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_util_functionspace_helpers_parallel
                  MPI     2
                  SOURCES test/util/FunctionSpaceHelpers.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_fullgmres
                  SOURCES test/assimilation/FullGMRES.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_rotmat
                  SOURCES test/assimilation/rotmat.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_solvematrixequation
                  SOURCES test/assimilation/SolveMatrixEquation.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_spectrallmp
                  SOURCES test/assimilation/SpectralLMP.cc
                  ARGS    "test/testinput/spectrallmp.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_testvector3d
                  SOURCES test/assimilation/TestVector3D.cc
                  ARGS    "test/testinput/empty.yaml"
                  LIBS    oops )

ecbuild_add_test( TARGET  test_assimilation_tridiagsolve
                  SOURCES test/assimilation/TriDiagSolve.cc
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

add_fctest( TARGET  test_util_missing_values_f.x
            SOURCES test/util/missing_values.F90
            LINKER_LANGUAGE Fortran
            CONDITION HAVE_FCTEST
            LIBS    oops )

endif()
