/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ModelTrajectory.h"

#include "eckit/exception/Exceptions.h"
#include "lorenz95/FieldL95.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ModelTrajectory::ModelTrajectory(const bool ltraj) : ltraj_(ltraj), traj_() {}
// -----------------------------------------------------------------------------
ModelTrajectory::~ModelTrajectory() {}
// -----------------------------------------------------------------------------
void ModelTrajectory::set(const FieldL95 & xx) {
  if (ltraj_) traj_.push_back(new FieldL95(xx));
}
// -----------------------------------------------------------------------------
const FieldL95 & ModelTrajectory::get(const int ii) const {
  ASSERT(ltraj_);
  ASSERT(traj_.size() == 4);
  ASSERT(1 <= ii && ii <= 4);
  return traj_[ii-1];
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95

