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
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/gatherPrint.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief State class used in oops; subclass of interface class interface::State.
///
/// \details
/// Handles additional MPI communicator parameter \p commTime_ in the constructors
/// (for MPI distribution in time, used in oops for 4DEnVar and weak-constraint 4DVar).
/// Adds communication through time to the following State methods:
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

  /// Constructor for specified \p resol and configuration \p conf
  ///
  /// \param conf
  ///   Configuration of either a single 3D state or a set of 3D states with different validity
  ///   times, to be used by different members of the MPI communicator in time
  State(const Geometry_ & resol, const eckit::Configuration & conf);

  /// Copies \p other State, changing its resolution to \p geometry
  State(const Geometry_ & resol, const State & other);

  /// Copies \p other State, changing its variables to \p vars
  State(const Variables & vars, const State & other);

  State(const State &);
  State & operator=(const State &);

  /// Accessor to geometry associated with this State
  const Geometry_ & geometry() const {return resol_;}

  /// Accessors to the FieldSet3D
  const FieldSet3D & fieldSet() const;
  FieldSet3D & fieldSet();
  void synchronizeFields();

  /// Write
  void write(const eckit::Configuration &) const;

 private:
  const Geometry_ & resol_;
};

// =============================================================================

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const Variables & vars,
                    const util::DateTime & time) :
  interface::State<MODEL>(resol, vars, time), resol_(resol)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const eckit::Configuration & conf) :
  interface::State<MODEL>(resol, conf), resol_(resol)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const State & other) :
  interface::State<MODEL>(resol, other), resol_(resol)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Variables & vars, const State & other) :
  interface::State<MODEL>(vars, other), resol_(other.resol_)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const State & other) :
  interface::State<MODEL>(other), resol_(other.resol_)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> & State<MODEL>::operator=(const State & rhs) {
  ASSERT(resol_ == rhs.resol_);
  interface::State<MODEL>::operator=(rhs);
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const FieldSet3D & State<MODEL>::fieldSet() const {
  if (!interface::State<MODEL>::fset_) {
    interface::State<MODEL>::fset_.reset(new FieldSet3D(this->validTime(), resol_.getComm()));
  }
  if (interface::State<MODEL>::fset_->empty()) {
    this->toFieldSet(interface::State<MODEL>::fset_->fieldSet());
    for (const auto & field : *interface::State<MODEL>::fset_) {
      ASSERT_MSG(field.rank() == 2,
                 "OOPS expects the model's State::toFieldSet method to return rank-2 fields,"
                 " but field " + field.name() + " has rank = " + std::to_string(field.rank()));
    }
  }
  return *interface::State<MODEL>::fset_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
FieldSet3D & State<MODEL>::fieldSet() {
  if (!interface::State<MODEL>::fset_) {
    interface::State<MODEL>::fset_.reset(new FieldSet3D(this->validTime(), resol_.getComm()));
  }
  if (interface::State<MODEL>::fset_->empty()) {
    this->toFieldSet(interface::State<MODEL>::fset_->fieldSet());
    for (const auto & field : *interface::State<MODEL>::fset_) {
      ASSERT_MSG(field.rank() == 2,
                 "OOPS expects the model's State::toFieldSet method to return rank-2 fields,"
                 " but field " + field.name() + " has rank = " + std::to_string(field.rank()));
    }
  }
  return *interface::State<MODEL>::fset_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::synchronizeFields() {
  // TODO(JEDI core team): remove this method when accessors are fully implemented
  ASSERT(interface::State<MODEL>::fset_);
  ASSERT(!interface::State<MODEL>::fset_->empty());
  this->fromFieldSet(interface::State<MODEL>::fset_->fieldSet());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::write(const eckit::Configuration & conf) const {
  const bool dateCols = conf.getBool("date colons", true);

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
      if (dateCols) {
        prefix += "." + antime.toString();
      } else {
        prefix += "." + antime.toStringIO();
      }
      const util::Duration step = this->validTime() - antime;
      prefix += "." + step.toString();
    }

    if (type == "an") {
      if (dateCols) {
        prefix += "." + this->validTime().toString();
      } else {
        prefix += "." + this->validTime().toStringIO();
      }
    }

    eckit::LocalConfiguration preconf(conf);
    if (conf.has("date colons")) preconf.set("date colons", dateCols);
    preconf.set("prefix", prefix);
    interface::State<MODEL>::write(preconf);

  } else {
    interface::State<MODEL>::write(conf);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATE_H_
