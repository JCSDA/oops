/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/DataSetBase.h"
#include "oops/base/Increment4D.h"
#include "oops/base/State4D.h"

namespace oops {

// -----------------------------------------------------------------------------
class FieldSet4D : public DataSetBase<FieldSet3D, atlas::FunctionSpace> {
  typedef DataSetBase<FieldSet3D, atlas::FunctionSpace> Base_;
 public:
  template<typename MODEL> FieldSet4D(const State4D<MODEL> &);
  template<typename MODEL> FieldSet4D(const Increment4D<MODEL> &);
 private:
  std::string classname() const {return "FieldSet4D";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSet4D::FieldSet4D(const State4D<MODEL> & state4d)
  : Base_(state4d.times(), state4d.commTime(), {0}, oops::mpi::myself()) {
  for (size_t jj = 0; jj < state4d.size(); ++jj) {
    this->dataset().emplace_back(new FieldSet3D(state4d[jj].fieldSet(),
                                                state4d[jj].validTime(),
                                                state4d.geometry().getComm()));
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSet4D::FieldSet4D(const Increment4D<MODEL> & inc4d)
  : Base_(inc4d.times(), inc4d.commTime(), {0}, oops::mpi::myself()) {
  for (size_t jj = 0; jj < inc4d.size(); ++jj) {
    this->dataset().emplace_back(new FieldSet3D(inc4d[jj].fieldSet(),
                                                inc4d[jj].validTime(),
                                                inc4d.geometry().getComm()));
  }
}

}  // namespace oops
