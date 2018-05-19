/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_MODELTRAJECTORY_H_
#define LORENZ95_MODELTRAJECTORY_H_

#include <string>

#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/util/ObjectCounter.h"

namespace lorenz95 {
  class FieldL95;

/// L95 model trajectory

// -----------------------------------------------------------------------------
class ModelTrajectory: private util::ObjectCounter<ModelTrajectory> {
 public:
  static const std::string classname() {return "lorenz95::ModelTrajectory";}

/// Constructor, destructor
  explicit ModelTrajectory(const bool ltraj = true);
  ~ModelTrajectory();

/// Save trajectory
  void set(const FieldL95 &);

/// Get trajectory
  const FieldL95 & get(const int) const;

 private:
  const bool ltraj_;
  boost::ptr_vector<FieldL95> traj_;
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_MODELTRAJECTORY_H_
