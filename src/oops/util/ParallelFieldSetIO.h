/*
 * (C) Crown copyright 2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARALLELFIELDSETIO_H_
#define OOPS_UTIL_PARALLELFIELDSETIO_H_

#include <algorithm>
#include <string>
#include <vector>

#include "atlas/redistribution/Redistribution.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/AtlasArrayUtil.h"

namespace util {

// -------------------------------------------------------------------------------------------------

class ParallelFieldSetIO {
 public:
    static const std::string classname() { return "oops::util::ParallelFieldSetIO"; }

    enum class Mode { Read, Write, ReadWrite };

    explicit ParallelFieldSetIO(const atlas::FunctionSpace&,
                                const std::string&,  // TODO(tom-j-h): remove after have Atlas #264
                                const Mode mode = Mode::ReadWrite);

    void write(const atlas::FieldSet&, const std::string&) const;
    void write(const atlas::Field&, const std::string&) const;

    void read(atlas::FieldSet&, const std::string&) const;
    void read(atlas::Field&, const std::string&) const;

 private:
    atlas::FunctionSpace setupFunctionSpace(const atlas::FunctionSpace&, const std::string&) const;

    atlas::FieldSet ioFieldSet(const atlas::FieldSet&) const;

    template<int Rank>
    std::vector<size_t> starts() const;

    template<int Rank>
    std::vector<size_t> counts(const atlas::Field&) const;

    void writeFieldByTypeAndRank(const atlas::Field&, std::vector<int>&, int&) const;

    void readFieldByTypeAndRank(atlas::Field&, const std::vector<int>&, const int) const;

    template <typename Value, int Rank>
    void writeField(const atlas::Field&, std::vector<int>&, int&) const;

    template <typename Value, int Rank>
    void readField(atlas::Field&, const std::vector<int>&, const int) const;

    template <typename Value>
    void dispatchWriteField(const atlas::Field&, std::vector<int>&, int&) const;

    template <typename Value>
    void dispatchReadField(atlas::Field&, const std::vector<int>&, const int) const;

    mutable atlas::idx_t gridSize_;  // TODO(tom-j-h): remove once Atlas #264 available, replace
                                     // all uses with functionSpace_.grid().size()
    const atlas::FunctionSpace functionSpace_;
    const atlas::Redistribution redistributeWrite_;
    const atlas::Redistribution redistributeRead_;
    const size_t startGidx_;
    const size_t countGidx_;
};

// -------------------------------------------------------------------------------------------------

template <int Rank>
std::vector<size_t> ParallelFieldSetIO::starts() const {
    std::vector<size_t> starts(Rank, 0);
    starts[0] = startGidx_;  // order data in file by global index
    return starts;
}

// -------------------------------------------------------------------------------------------------

template <int Rank>
std::vector<size_t> ParallelFieldSetIO::counts(const atlas::Field& field) const {
    std::vector<size_t> counts(Rank);
    const std::vector<atlas::idx_t> shape(field.shape());
    std::transform(shape.begin(), shape.end(), counts.begin(),
                   [](auto x){ return static_cast<size_t>(x); });  // netCDF API requires size_t
    counts[0] = countGidx_;  // ignore ghost points
    return counts;
}

// -------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void ParallelFieldSetIO::writeField(const atlas::Field& field,
                                    std::vector<int>& netcdfGeneralIDs,
                                    int& netcdfVarID) const {
    auto view = atlas::array::make_view<const Value, Rank>(field);
    atlasArrayWriteData(netcdfGeneralIDs, netcdfVarID, starts<Rank>(), counts<Rank>(field), view);
}

// -------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void ParallelFieldSetIO::readField(atlas::Field& field,
                                   const std::vector<int>& netcdfGeneralIDs,
                                   const int netcdfVarID) const {
    auto view = atlas::array::make_view<Value, Rank>(field);
    atlasArrayReadData(netcdfGeneralIDs, netcdfVarID, starts<Rank>(), counts<Rank>(field), view);
}

// -------------------------------------------------------------------------------------------------

template <typename Value>
void ParallelFieldSetIO::dispatchWriteField(const atlas::Field& field,
                                            std::vector<int>& netcdfGeneralIDs,
                                            int& netcdfVarID) const {
    switch (field.rank()) {
        case 1:
            writeField<Value, 1>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case 2:
            writeField<Value, 2>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case 3:
            writeField<Value, 3>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        default:
            ABORT("util::ParallelFieldSetIO::dispatchWriteField: Rank of Field must be 3 or less");
    }
}

// -------------------------------------------------------------------------------------------------

template <typename Value>
void ParallelFieldSetIO::dispatchReadField(atlas::Field& field,
                                           const std::vector<int>& netcdfGeneralIDs,
                                           const int netcdfVarID) const {
    switch (field.rank()) {
        case 1:
            readField<Value, 1>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case 2:
            readField<Value, 2>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        case 3:
            readField<Value, 3>(field, netcdfGeneralIDs, netcdfVarID);
            break;
        default:
            ABORT("util::ParallelFieldSetIO::dispatchReadField: Rank of Field must be 3 or less");
    }
}

// -------------------------------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_PARALLELFIELDSETIO_H_
