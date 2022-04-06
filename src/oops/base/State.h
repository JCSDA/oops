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

#include "oops/base/StateParametersND.h"
#include "oops/interface/State.h"
#include "oops/util/gatherPrint.h"

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
  typedef typename interface::State<MODEL>::Parameters_ Parameters_;
  typedef typename interface::State<MODEL>::WriteParameters_ WriteParameters_;

  /// Configuration options of either a single 3D model state  or a set of 3D states valid at
  /// different times, each used by a different member of the MPI communicator in time.
  typedef StateParametersND<MODEL> ParametersND_;

  /// Constructor for specified \p resol, with \p vars, valid at \p time
  State(const Geometry_ & resol, const Variables & vars, const util::DateTime & time);

  /// Constructor for specified \p resol and parameters \p params
  ///
  /// \param params
  ///   Parameters of either a single 3D state or a set of 3D states with different validity times,
  ///   to be used by different members of the MPI communicator in time
  State(const Geometry_ & resol, const ParametersND_ & params);

  /// Constructor for specified \p resol and configuration \p conf
  ///
  /// \param conf
  ///   Configuration of either a single 3D state or a set of 3D states with different validity
  ///   times, to be used by different members of the MPI communicator in time
  State(const Geometry_ & resol, const eckit::Configuration & conf);

  /// Constructor for specified \p resol and parameters \p params
  ///
  /// \param params
  ///   Parameters of a single 3D state
  State(const Geometry_ & resol, const Parameters_ & params);
  /// Copies \p other State, changing its resolution to \p geometry
  State(const Geometry_ & resol, const State & other);

  State(const State &);
  State & operator=(const State &);

  /// Accessor to geometry associated with this State
  const Geometry_ & geometry() const {return resol_;}

  // Get values as Atlas FieldSet
  const atlas::FieldSet & fieldSet() const;

  /// Norm (used in tests)
  double norm() const;

  /// Write
  void write(const eckit::Configuration &) const;
  void write(const WriteParameters_ &) const;

 private:
  const Geometry_ & resol_;
  const eckit::mpi::Comm * commTime_;  /// pointer to the MPI communicator in time

  void print(std::ostream &) const override;
};

// =============================================================================

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const Variables & vars,
                    const util::DateTime & time) :
  interface::State<MODEL>(resol, vars, time), resol_(resol), commTime_(&resol.timeComm())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const ParametersND_ & paramsND) :
  // The call to at() will return the parameters of the 3D state to be used by the current MPI rank
  State(resol, paramsND.at(resol.timeComm()))
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const Parameters_ & params) :
  interface::State<MODEL>(resol, params), resol_(resol), commTime_(&resol.timeComm())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const eckit::Configuration & conf) :
  State(resol, validateAndDeserialize<ParametersND_>(conf))
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const State & other) :
  interface::State<MODEL>(resol, other), resol_(resol), commTime_(&resol.timeComm())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const State & other) :
  interface::State<MODEL>(other.resol_, other), resol_(other.resol_), commTime_(&resol_.timeComm())
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> & State<MODEL>::operator=(const State & rhs) {
  ASSERT(resol_ == rhs.resol_);
  ASSERT(commTime_ == rhs.commTime_);
  interface::State<MODEL>::operator=(rhs);
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const atlas::FieldSet & State<MODEL>::fieldSet() const {
  if (!interface::State<MODEL>::fset_) {
    interface::State<MODEL>::fset_ = std::make_unique<atlas::FieldSet>();
    this->getFieldSet(this->variables(), *interface::State<MODEL>::fset_);
  }
  return *interface::State<MODEL>::fset_;
}

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
void State<MODEL>::write(const eckit::Configuration & conf) const {
  if (conf.has("type") && conf.has("exp") && !conf.has("prefix")) {
    const std::string type = conf.getString("type");
    std::string prefix = conf.getString("exp") + "." + type;

    if (type == "ens") {
      if (!conf.has("member"))
        throw eckit::BadValue("'member' was not set in the parameters passed to write() "
                              "even though 'type' was set to '" + type + "'", Here());
      prefix += "." + std::to_string(conf.getInt("member"));
    }

    if (type == "fc" || type == "ens") {
      if (!conf.has("date"))
        throw eckit::BadValue("'date' was not set in the parameters passed to write() "
                              "even though 'type' was set to '" + type + "'", Here());
      const util::DateTime antime(conf.getString("date"));
      prefix += "." + antime.toString();
      const util::Duration step = this->validTime() - antime;
      prefix += "." + step.toString();
    }

    if (type == "an") {
      prefix += "." + this->validTime().toString();
    }

    eckit::LocalConfiguration preconf(conf);
    preconf.set("prefix", prefix);
    interface::State<MODEL>::write(preconf);
  } else {
    interface::State<MODEL>::write(conf);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::write(const WriteParameters_ & params) const {
  eckit::LocalConfiguration conf = params.toConfiguration();
  this->write(conf);
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
