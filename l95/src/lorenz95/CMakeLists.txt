set( _l95_srcs
        AnalyticInit.cc
        AnalyticInit.h
        BackgroundCheck.cc
        BackgroundCheck.h
        ErrorCovarianceL95.cc
        ErrorCovarianceL95.h
        FieldL95.cc
        FieldL95.h
        GetValuesL95.cc
        GetValuesL95.h
        GetValuesTLAD.cc
        GetValuesTLAD.h
        GomL95.cc
        GomL95.h
        GomL95.cc
        GomL95.h
        IncrementL95.cc
        IncrementL95.h
        instantiateL95ChangeVarFactory.h
        instantiateLocalizationFactory.h
        Iterator.cc
        Iterator.h
        L95Traits.h
        LocsL95.cc
        LocsL95.h
        ModelBias.cc
        ModelBias.h
        ModelBiasCorrection.cc
        ModelBiasCorrection.h
        ModelBiasCovariance.cc
        ModelBiasCovariance.h
        ModelL95.h
        ModelL95.cc
        ModelTrajectory.cc
        ModelTrajectory.h
        ObsBias.cc
        ObsBias.h
        ObsBiasCorrection.cc
        ObsBiasCorrection.h
        ObsBiasCovariance.cc
        ObsBiasCovariance.h
        ObsBiasParameters.h
        ObsDiags1D.h
        ObsIterator.cc
        ObsIterator.h
        ObsLocGC99.cc
        ObsLocGC99.h
        ObsTable.cc
        ObsTable.h
        ObsVec1D.cc
        ObsVec1D.h
        ObsData1D.h
        ObservationL95.cc
        ObservationL95.h
        ObservationTLAD.cc
        ObservationTLAD.h
        QCmanager.cc
        QCmanager.h
        Resolution.cc
        Resolution.h
        StateL95.cc
        StateL95.h
        TLML95.cc
        TLML95.h
        LocalizationMatrixL95.cc
        LocalizationMatrixL95.h )

ecbuild_add_library( TARGET lorenz95
                     SOURCES ${_l95_srcs}
                     PUBLIC_LIBS oops
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/l95
                     LINKER_LANGUAGE CXX )

#TODO: Reenable cxx_std_14 when supported by intel jedi-stack modules (jedi-stack #125)
target_compile_features( lorenz95 PUBLIC cxx_std_11 )
target_include_directories( lorenz95 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> )
