/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/base/DataSetBase.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/IncrementSet.h"

namespace oops {

// -----------------------------------------------------------------------------
class FieldSets : public DataSetBase<FieldSet3D, atlas::FunctionSpace> {
  typedef DataSetBase<FieldSet3D, atlas::FunctionSpace> Base_;

 public:
  /// @brief Creates a FieldSets from the IncrementSet. On creation fieldsets are
  ///        shared between Increment::fieldSet() and fieldsets in FieldSet4D.
  template<typename MODEL> FieldSets(const IncrementSet<MODEL> &);

  FieldSets(const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
            const std::vector<int> &, const eckit::mpi::Comm &);

  FieldSets(const atlas::FunctionSpace &,
            const Variables &,
            const std::vector<util::DateTime> &,
            const eckit::Configuration &,
            const eckit::mpi::Comm &,
            const eckit::mpi::Comm & = oops::mpi::myself(),
            const eckit::mpi::Comm & = oops::mpi::myself());

  /// @brief  Multiplies each FieldSet3D in this FieldSets with the \p other.
  FieldSets & operator*=(const oops::FieldSet3D & other);
  FieldSets & operator*=(const double zz);

 private:
  std::string classname() const {return "FieldSets";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSets::FieldSets(const IncrementSet<MODEL> & inc4dens)
  : Base_(inc4dens.times(), inc4dens.commTime(), inc4dens.members(), inc4dens.commEns()) {
  for (size_t jj = 0; jj < inc4dens.size(); ++jj) {
    this->dataset().emplace_back(std::make_unique<FieldSet3D>(inc4dens[jj].validTime(),
                                                inc4dens.geometry().getComm()));
    this->dataset()[jj]->shallowCopy(inc4dens[jj].fieldSet());
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
