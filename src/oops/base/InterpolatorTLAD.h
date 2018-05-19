/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_INTERPOLATORTLAD_H_
#define OOPS_BASE_INTERPOLATORTLAD_H_

#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/interface/InterpolatorTraj.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class InterpolatorTLAD {
  typedef InterpolatorTraj<MODEL>    InterpolatorTraj_;

 public:
  explicit InterpolatorTLAD(const std::size_t nn): traj_() {
    for (size_t jj = 0; jj < nn; ++jj) {
      boost::shared_ptr<InterpolatorTraj_> pp(new InterpolatorTraj_());
      traj_.push_back(pp);
    }
    Log::trace() << "InterpolatorTLAD::InterpolatorTLAD size = " << nn << std::endl;
  }
  ~InterpolatorTLAD() {}

  InterpolatorTraj_ & operator[](const std::size_t ii) {return *traj_.at(ii);}
  const InterpolatorTraj_ & operator[](const std::size_t ii) const {return *traj_.at(ii);}

 private:
  std::vector<boost::shared_ptr<InterpolatorTraj_> > traj_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INTERPOLATORTLAD_H_
