list( APPEND oops_qg_src_files
AnalyticInit.cc
AnalyticInit.h
ChangeVarQG.cc
ChangeVarQG.h
ChangeVarTLADQG.cc
ChangeVarTLADQG.h
ErrorCovarianceQG.cc
ErrorCovarianceQG.h
FieldsQG.cc
FieldsQG.h
GeometryQG.cc
GeometryQG.h
GeometryQGIterator.cc
GeometryQGIterator.h
GomQG.cc
GomQG.h
IncrementQG.cc
IncrementQG.h
InterpolatorQG.cc
InterpolatorQG.h
LocationsQG.cc
LocationsQG.h
LocalizationMatrixQG.cc
LocalizationMatrixQG.h
ModelBias.h
ModelBiasCovariance.h
ModelBiasIncrement.h
ModelData.cc
ModelData.h
ModelQG.cc
ModelQG.h
ObsBias.cc
ObsBias.h
ObsBiasCovariance.cc
ObsBiasCovariance.h
ObsBiasIncrement.cc
ObsBiasIncrement.h
ObsBiasPreconditioner.cc
ObsBiasPreconditioner.h
ObsDataQG.h
ObsDiagsQG.h
ObsFilter.cc
ObsFilter.h
ObsIteratorQG.cc
ObsIteratorQG.h
ObsLocQG.cc
ObsLocQG.h
ObsOpBaseQG.cc
ObsOpBaseQG.h
ObsOpBaseTLAD.cc
ObsOpBaseTLAD.h
ObsOperatorQG.cc
ObsOperatorQG.h
ObsOperatorTLAD.cc
ObsOperatorTLAD.h
ObsSpaceQG.cc
ObsSpaceQG.h
ObsStreamQG.cc
ObsStreamQG.h
ObsStreamTLAD.cc
ObsStreamTLAD.h
ObsVecQG.cc
ObsVecQG.h
ObsWSpeedQG.cc
ObsWSpeedQG.h
ObsWSpeedTLAD.cc
ObsWSpeedTLAD.h
ObsWindQG.cc
ObsWindQG.h
ObsWindTLAD.cc
ObsWindTLAD.h
QgFortran.h
QgTraits.h
StateQG.cc
StateQG.h
TlmQG.cc
TlmQG.h
fft_f.cc
fft_f.h
fft_mod.F90
instantiateQgNormFactory.h
instantiateQgLocalizationFactory.h
qg_advect_q_mod.F90
qg_change_var_interface.F90
qg_change_var_mod.F90
qg_constants_mod.F90
qg_convert_q_to_x_mod.F90
qg_convert_x_to_q_mod.F90
qg_convert_x_to_u_mod.F90
qg_convert_x_to_v_mod.F90
qg_differential_solver_mod.F90
qg_error_covariance_interface.F90
qg_error_covariance_mod.F90
qg_fields_interface.F90
qg_fields_mod.F90
qg_geom_interface.F90
qg_geom_mod.F90
qg_geom_iter_interface.F90
qg_geom_iter_mod.F90
qg_getvalues_interp_interface.F90
qg_getvalues_interp_mod.F90
qg_gom_interface.F90
qg_gom_mod.F90
qg_interp_mod.F90
qg_locs_f.cc
qg_locs_f.h
qg_locs_mod.F90
qg_model_interface.F90
qg_model_mod.F90
qg_obsdb_interface.F90
qg_obsdb_mod.F90
qg_obsvec_interface.F90
qg_obsvec_mod.F90
qg_obsvec_random_f.h
qg_obsvec_random_f.cc
qg_projection_mod.F90
qg_stream_interface.F90
qg_stream_mod.F90
qg_tools_mod.F90
qg_wind_interface.F90
qg_wind_mod.F90
qg_wspeed_interface.F90
qg_wspeed_mod.F90
)

ecbuild_add_library( TARGET   qg
                     SOURCES  ${oops_qg_src_files}
                     PUBLIC_LIBS     oops
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/qg
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

#Configure include directory layout for build-tree to match install-tree
set(QG_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(qg_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${QG_BUILD_DIR_INCLUDE_PATH}/oops"
                                    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${QG_BUILD_DIR_INCLUDE_PATH}/oops/qg)
add_dependencies(qg qg_headers)
target_include_directories(qg PUBLIC $<BUILD_INTERFACE:${QG_BUILD_DIR_INCLUDE_PATH}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../..>)

#Fortran module output directory for build interface
set(QG_MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION}/qg)
set_target_properties(qg PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${QG_MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${QG_MODULE_DIR}/ DESTINATION ${QG_MODULE_DIR})
target_include_directories(qg INTERFACE
                                $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${QG_MODULE_DIR}>
                                $<INSTALL_INTERFACE:${QG_MODULE_DIR}>)
