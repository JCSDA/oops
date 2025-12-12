/*
 * (C) Crown copyright 2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/ParallelFieldSetIO.h"

#include <string>
#include <vector>

#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Timer.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace util {

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace ParallelFieldSetIO::setupFunctionSpace(
    const atlas::FunctionSpace& nativeFunctionSpace, const std::string& gridName) const {
    atlas::Grid grid{};
    atlas::grid::Partitioner partitioner{};
    atlas::Mesh mesh{};
    atlas::FunctionSpace functionSpace{};
    atlas::FieldSet fieldSet{};

    eckit::LocalConfiguration functionSpaceConfig;
    // TODO(tom-j-h): once Atlas #264 available, call nativeFunctionSpace.grid().name() instead,
    // and remove gridName argument from this function and class constructor
    functionSpaceConfig.set("grid.name", gridName);
    functionSpaceConfig.set("function space", nativeFunctionSpace.type());
    functionSpaceConfig.set("halo", 0);
    // This partitioner splits up the Grid ~evenly across all PEs, and orders points by global_index
    functionSpaceConfig.set("partitioner", "equal_bands");

    util::setupFunctionSpace(atlas::mpi::comm(), functionSpaceConfig, grid, partitioner, mesh,
                             functionSpace, fieldSet);

    gridSize_ = grid.size();  // TODO(tom-j-h): remove once Atlas #264 available

    return functionSpace;
}

// -------------------------------------------------------------------------------------------------

ParallelFieldSetIO::ParallelFieldSetIO(const atlas::FunctionSpace& nativeFunctionSpace,
                                       const std::string& gridName,
                                       const Mode mode)
    : functionSpace_(setupFunctionSpace(nativeFunctionSpace, gridName)),
      redistributeWrite_((mode == Mode::Write || mode == Mode::ReadWrite) ?
          atlas::Redistribution(nativeFunctionSpace, functionSpace_) : atlas::Redistribution()),
      redistributeRead_((mode == Mode::Read || mode == Mode::ReadWrite) ?
          atlas::Redistribution(functionSpace_, nativeFunctionSpace) : atlas::Redistribution())
    {}

// -------------------------------------------------------------------------------------------------

atlas::FieldSet ParallelFieldSetIO::ioFieldSet(const atlas::FieldSet& nativeFieldSet) const {
    atlas::FieldSet ioFieldSet;
    for (const auto& field : nativeFieldSet) {
        const std::vector<atlas::idx_t> shape(field.shape());
        ASSERT_MSG(shape.size() < 4,
                   "util::ParallelFieldSetIO::ioFieldSet: Rank is greater than 3");
        // TODO(tom-j-h): Add following check once code in Atlas #264 is available
        // ASSERT_MSG(field.functionspace().grid() == functionSpace_.grid(),
        //            "util::ParallelFieldSetIO::ioFieldSet: Field's Grid doesn't match I/O Grid");
        auto fieldSetConfig = atlas::option::name(field.name());
        fieldSetConfig.set(atlas::option::halo(0));
        if (shape.size() > 1) fieldSetConfig.set(atlas::option::levels(shape[1]));
        if (shape.size() > 2) fieldSetConfig.set(atlas::option::vector(shape[2]));
        switch (field.datatype().kind()) {
            case atlas::array::DataType::KIND_INT32:
                ioFieldSet.add(functionSpace_.createField<int>(fieldSetConfig));
                break;
            case atlas::array::DataType::KIND_REAL32:
                // Note - although this class can deal with 32-bit floating point data,
                // AtlasArrayUtil cannot yet, and will throw an exception.
                // TODO(tom-j-h): Add further datatypes to AtlasArrayUtil.
                ioFieldSet.add(functionSpace_.createField<float>(fieldSetConfig));
                break;
            case atlas::array::DataType::KIND_REAL64:
                ioFieldSet.add(functionSpace_.createField<double>(fieldSetConfig));
                break;
            default:
                ABORT("util::ParallelFieldSetIO::ioFieldSet: invalid datatype for Atlas Field");
        }
    }
    return ioFieldSet;
}

// -------------------------------------------------------------------------------------------------

void ParallelFieldSetIO::writeFieldByTypeAndRank(const atlas::Field& field,
                                                 std::vector<int>& netcdfGeneralIDs,
                                                 int& netcdfVarID) const {
    switch (field.datatype().kind()) {
        case atlas::array::DataType::KIND_INT32:
            dispatchWriteField<int>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case atlas::array::DataType::KIND_REAL32:
            dispatchWriteField<float>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case atlas::array::DataType::KIND_REAL64:
            dispatchWriteField<double>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        default:
            ABORT("util::ParallelFieldSetIO::writeFieldByTypeAndRank: "
                  "invalid datatype for Atlas Field");
    }
}

// -------------------------------------------------------------------------------------------------

void ParallelFieldSetIO::readFieldByTypeAndRank(atlas::Field& field,
                                                const std::vector<int>& netcdfGeneralIDs,
                                                const int netcdfVarID) const {
    switch (field.datatype().kind()) {
        case atlas::array::DataType::KIND_INT32:
            dispatchReadField<int>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case atlas::array::DataType::KIND_REAL32:
            dispatchReadField<float>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case atlas::array::DataType::KIND_REAL64:
            dispatchReadField<double>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        default:
            ABORT("util::ParallelFieldSetIO::readFieldByTypeAndRank: "
                  "invalid datatype for Atlas Field");
    }
}

// -------------------------------------------------------------------------------------------------

void ParallelFieldSetIO::write(const atlas::FieldSet& source, const std::string& ncfilepath) const {
    ASSERT_MSG(redistributeWrite_,
               "util::ParallelFieldSetIO::write: redistributeWrite_ not initialised");
    util::Timer timer(classname(), "write");

    // Redistribute data to Fields where it is ordered by global_index
    atlas::FieldSet output = ioFieldSet(source);
    redistributeWrite_.execute(source, output);

    // Obtain dimension information to pass to I/O routines
    // TODO(tom-j-h): once Atlas #264 available, call functionSpace_.grid().size() rather than using
    // gridSize_
    std::vector<atlas::idx_t> dimSizes = {gridSize_}, r2DimSizes, r3DimSizes;
    std::vector<std::string> dimNames = {"global_index"}, r2DimNames, r3DimNames;
    std::vector<std::vector<std::string>> dimNamesForEveryVar;
    dimNamesForEveryVar.reserve(source.size());
    size_t r2DimCounter = 0, r3DimCounter = 0;
    for (const auto& field : source) {
        const std::vector<atlas::idx_t> shape(field.shape());
        std::vector<std::string> dimNamesForThisVar(shape.size());
        dimNamesForThisVar[0] = dimNames[0];  // first dimension of each Field is always Grid size
        if (shape.size() > 1) {
            const bool newNumberOfLevels = std::find(
                r2DimSizes.begin(), r2DimSizes.end(), shape[1]) == r2DimSizes.end();
            if (newNumberOfLevels) {
                // If no previous Field has this number of levels, make a new dimension
                r2DimSizes.push_back(shape[1]);
                const std::string dimName = "r2d" + std::to_string(r2DimCounter);
                r2DimNames.push_back(dimName);
                dimNamesForThisVar[1] = dimName;
                r2DimCounter++;
            } else {
                // Otherwise use the existing dimension
                const auto index = std::distance(
                    r2DimSizes.begin(), std::find(r2DimSizes.begin(), r2DimSizes.end(), shape[1]));
                dimNamesForThisVar[1] = r2DimNames[index];
            }
        }
        if (shape.size() > 2) {
            const bool newVectorLength = std::find(
                r3DimSizes.begin(), r3DimSizes.end(), shape[2]) == r3DimSizes.end();
            if (newVectorLength) {
                // If no previous Field has this 3rd rank size, make a new dimension
                r3DimSizes.push_back(shape[2]);
                const std::string dimName = "r3d" + std::to_string(r3DimCounter);
                r3DimNames.push_back(dimName);
                dimNamesForThisVar[2] = dimName;
                r3DimCounter++;
            } else {
                // Otherwise use the existing dimension
                const auto index = std::distance(
                    r3DimSizes.begin(), std::find(r3DimSizes.begin(), r3DimSizes.end(), shape[2]));
                dimNamesForThisVar[2] = r3DimNames[index];
            }
        }
        dimNamesForEveryVar.push_back(dimNamesForThisVar);
    }
    dimSizes.insert(dimSizes.end(), r2DimSizes.begin(), r2DimSizes.end());
    dimSizes.insert(dimSizes.end(), r3DimSizes.begin(), r3DimSizes.end());
    dimNames.insert(dimNames.end(), r2DimNames.begin(), r2DimNames.end());
    dimNames.insert(dimNames.end(), r3DimNames.begin(), r3DimNames.end());

    // Write header
    std::vector<oops::Variable> variablesVec;
    variablesVec.reserve(output.size());
    for (auto f = 0; f < output.size(); f++) {
        oops::ModelDataType dt;
        switch (output[f].datatype().kind()) {
            case atlas::array::DataType::KIND_INT32:
                dt = oops::ModelDataType::Int32;
                break;
            case atlas::array::DataType::KIND_REAL32:
                dt = oops::ModelDataType::Real32;
                break;
            case atlas::array::DataType::KIND_REAL64:
                dt = oops::ModelDataType::Real64;
                break;
            default:
                ABORT("util::ParallelFieldSetIO::write: invalid datatype for Atlas Field");
        }
        variablesVec.push_back(oops::Variable(
            output[f].name(),
            oops::VariableMetaData(oops::defaultVerticalStagger, dt, oops::defaultVariableDomain),
            output[f].levels()));
    }
    const auto variables = oops::Variables(variablesVec);
    std::vector<int> netcdfGeneralIDs, netcdfDimIDs, netcdfVarIDs;
    std::vector<std::vector<int>> netcdfDimVarIDs;
    const eckit::LocalConfiguration config;
    oops::Log::info() << "util::ParallelFieldSetIO::write: writing file "
                      << ncfilepath << std::endl;
    atlasArrayWriteHeader(ncfilepath, dimNames, dimSizes, variables, dimNamesForEveryVar,
                          config, netcdfGeneralIDs, netcdfDimIDs, netcdfVarIDs, netcdfDimVarIDs,
                          true);

    // Write Fields
    for (auto f = 0; f < output.size(); f++) {
        const auto& field = output[f];
        writeFieldByTypeAndRank(field, netcdfGeneralIDs, netcdfVarIDs[f]);
    }

    // Close file
    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ERR(retval);
    oops::Log::info() << "util::ParallelFieldSetIO::write: finished writing file "
                      << ncfilepath << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ParallelFieldSetIO::write(const atlas::Field& source, const std::string& ncfilepath) const {
    atlas::FieldSet fieldSet;
    fieldSet.add(source);
    write(fieldSet, ncfilepath);
}

// -------------------------------------------------------------------------------------------------

void ParallelFieldSetIO::read(atlas::FieldSet& target, const std::string& ncfilepath) const {
    util::Timer timer(classname(), "read");

    ASSERT_MSG(redistributeRead_,
               "util::ParallelFieldSetIO::read: redistributeRead_ not initialised");

    // Set up Fields where data is ordered the same as in the file
    atlas::FieldSet input = ioFieldSet(target);

    // Set up variables to pass to I/O routines to obtain dimension information
    oops::Variables variables;
    std::vector<std::string> dimNames;
    std::vector<atlas::idx_t> dimSizes;
    std::vector<std::vector<std::string>> dimNamesForEveryVar;
    std::vector<int> netcdfGeneralIDs, netcdfDimIDs, netcdfVarIDs;
    std::vector<std::vector<int>> netcdfDimVarIDs;
    eckit::LocalConfiguration config;

    // Read header
    oops::Log::info() << "util::ParallelFieldSetIO::read: reading file "
                      << ncfilepath << std::endl;
    atlasArrayInquire(ncfilepath, dimNames, dimSizes, variables, dimNamesForEveryVar, config,
                      netcdfGeneralIDs, netcdfDimIDs, netcdfVarIDs, netcdfDimVarIDs, true);

    // Read Fields
    for (auto f = 0; f < input.size(); f++) {
        auto& field = input[f];
        const auto varID = variables.find(input.field_names()[f]);
          // in case we're not reading all variables
        readFieldByTypeAndRank(field, netcdfGeneralIDs, varID);
    }

    // Close file
    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ERR(retval);
    oops::Log::info() << "util::ParallelFieldSetIO::read: finished reading file "
                      << ncfilepath << std::endl;

    // Redistribute data to target Fields
    redistributeRead_.execute(input, target);
    target.set_dirty();
    target.haloExchange();
}

// -------------------------------------------------------------------------------------------------

void ParallelFieldSetIO::read(atlas::Field& target, const std::string& ncfilepath) const {
    atlas::FieldSet fieldSet;
    fieldSet.add(target);
    read(fieldSet, ncfilepath);
}

// -------------------------------------------------------------------------------------------------

}  // namespace util
