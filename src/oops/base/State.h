/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATE_H_
#define OOPS_BASE_STATE_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/interface/State.h"


namespace oops {

// -----------------------------------------------------------------------------
/// \brief State class used in oops; subclass of interface class interface::State.
///
/// \details
/// Handles additional MPI communicator parameter \p commTime_ in the constructors
/// (for MPI distribution in time, used in oops for 4DEnVar and weak-constraint 4DVar).
/// Adds communication through time to the following Increment methods:
/// - norm
/// - print

// -----------------------------------------------------------------------------
template <typename MODEL>
class State : public interface::State<MODEL> {
  typedef typename MODEL::State              State_;
  typedef Geometry<MODEL>                    Geometry_;

 public:
  /// Constructor for specified \p resol, with \p vars, valid at \p time
  State(const Geometry_ & resol, const Variables & vars, const util::DateTime & time);
  /// Constructor for specified \p resol and files read from \p conf
  State(const Geometry_ & resol, const eckit::Configuration & conf);
  /// Copies \p other State, changing its resolution to \p geometry
  State(const Geometry_ & resol, const State & other);

  /// Norm (used in tests)
  double norm() const;

 private:
  const eckit::mpi::Comm * commTime_;  /// pointer to the MPI communicator in time
  void print(std::ostream &) const override;
};

// =============================================================================

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const Variables & vars,
                    const util::DateTime & time) :
  interface::State<MODEL>(resol, vars, time), commTime_(&resol.timeComm())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const eckit::Configuration & conf) :
  interface::State<MODEL>(resol, conf), commTime_(&resol.timeComm())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const State & other) :
  interface::State<MODEL>(resol, other), commTime_(&resol.timeComm())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
double State<MODEL>::norm() const {
  double zz = interface::State<MODEL>::norm();
  zz *= zz;
  commTime_->allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  zz = sqrt(zz);
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::print(std::ostream & os) const {
  if (commTime_->size() > 1) {
    gatherPrint(os, this->state(), *commTime_);
  } else {
    os << this->state();
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATE_H_
