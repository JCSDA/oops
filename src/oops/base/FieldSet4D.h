/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/FieldSet3D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/Increment4D.h"
#include "oops/base/State4D.h"

namespace oops {

// -----------------------------------------------------------------------------
class FieldSet4D : public FieldSets {
 public:
  /// @brief Creates a FieldSet4D for specified times with specified time
  ///        communicator. Data are allocated in the ctor, all fieldsets are
  ///        empty.
  FieldSet4D(const std::vector<util::DateTime> & times,
             const eckit::mpi::Comm & commTime,
             const eckit::mpi::Comm & commGeom);
  /// @brief Creates a FieldSet4D with a single 3D fieldset in it.
  explicit FieldSet4D(const FieldSet3D &);
  /// @brief Creates a FieldSet4D from the State4D. On creation fieldsets are
  ///        shared between State::fieldSet() and fieldsets in FieldSet4D.
  template<typename MODEL> FieldSet4D(const State4D<MODEL> &);
  /// @brief Creates a FieldSet4D from the Increment4D. On creation fieldsets are
  ///        shared between Increment::fieldSet() and fieldsets in FieldSet4D.
  template<typename MODEL> FieldSet4D(const Increment4D<MODEL> &);

  /// @brief Initialize with a deep copy of \p iens ensemble member of \p other,
  ///        consistency check for valid time
  void deepCopy(const FieldSets & other, const size_t iens);

  void zero();
  FieldSet4D & operator+=(const FieldSet4D & other);
  FieldSet4D & operator*=(const FieldSet4D & other);
  /// @brief  Multiplies each FieldSet3D in this FieldSet4D with the \p other.
  FieldSet4D & operator*=(const FieldSet3D & other);
  FieldSet4D & operator*=(const double zz);
  /// @brief Computes dot product of this FieldSet4D with the \p other FieldSet4D
  ///        only for specified variables \p vars.
  double dot_product_with(const FieldSet4D & other, const Variables & vars) const;
  /// @brief Computes dot product of this FieldSet4D with the \p iens member
  ///        of the \p other FieldSets only for specified variables \p vars.
  double dot_product_with(const FieldSets & other, const size_t iens,
                          const Variables & vars) const;
  double norm() const;

 private:
  std::string classname() const {return "FieldSet4D";}
};

/// @brief Deep or shallow copy of the FieldSet4D. Using FieldSet4D copy ctor or assignment
///        operator results in deep copies.
FieldSet4D copyFieldSet4D(const FieldSet4D & other,
                          const bool & shallow = false);

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSet4D::FieldSet4D(const State4D<MODEL> & state4d)
  : FieldSets(state4d.times(), state4d.commTime(), {0}, oops::mpi::myself()) {
  for (size_t jj = 0; jj < state4d.size(); ++jj) {
    this->dataset().emplace_back(new FieldSet3D(state4d[jj].validTime(),
                                                state4d.geometry().getComm()));
    this->dataset()[jj]->shallowCopy(state4d[jj].fieldSet());
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSet4D::FieldSet4D(const Increment4D<MODEL> & inc4d)
  : FieldSets(inc4d.times(), inc4d.commTime(), {0}, oops::mpi::myself()) {
  for (size_t jj = 0; jj < inc4d.size(); ++jj) {
    this->dataset().emplace_back(new FieldSet3D(inc4d[jj].validTime(),
                                                inc4d.geometry().getComm()));
    this->dataset()[jj]->shallowCopy(inc4d[jj].fieldSet());
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
