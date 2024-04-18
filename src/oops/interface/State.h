/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_STATE_H_
#define OOPS_INTERFACE_STATE_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"
#include "oops/util/Timer.h"

namespace atlas {
  class FieldSet;
}

namespace oops {

namespace interface {

/// \brief Encapsulates the model state.
// -----------------------------------------------------------------------------

template <typename MODEL>
class State : public util::Printable,
              public util::Serializable,
              private util::ObjectCounter<State<MODEL> > {
  typedef typename MODEL::State            State_;
  typedef oops::Geometry<MODEL>            Geometry_;

 public:
  static const std::string classname() {return "oops::State";}

  /// Constructor for specified \p resol, with \p vars, valid at \p time
  State(const Geometry_ & resol, const Variables & vars, const util::DateTime & time);
  /// Constructor for specified \p resol and parameters \p params specifying e.g. a file to read
  /// or an analytic state to generate
  State(const Geometry_ & resol, const eckit::Configuration & conf);
  /// Copies \p other State, changing its resolution to \p geometry
  State(const Geometry_ & resol, const State & other);
  /// Copies \p other State, changing its variables to \p vars
  State(const Variables & vars, const State & other);
  /// Copy constructor
  State(const State &);
  /// Destructor (defined explicitly for timing and tracing)
  ~State();
  /// Assignment operator
  State & operator =(const State &);

  /// Accessor
  State_ & state() {if (fset_) {fset_->clear();} return *state_;}
  /// const accessor
  const State_ & state() const {return *state_;}

  /// Accessor to the time of this State
  const util::DateTime validTime() const {return state_->validTime();}
  /// Update this State's valid time by \p dt
  void updateTime(const util::Duration & dt) {state_->updateTime(dt);}

  /// Get run ID (used for generic PseudoModel class, set to -1 if not needed)
  size_t ID() const {return ID_;}
  /// Read this State from file
  void read(const eckit::Configuration &);
  /// Write this State out to file
  void write(const eckit::Configuration &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Accessor to variables associated with this State
  const Variables & variables() const;

  /// Zero out this State
  void zero();
  /// Accumulate (add \p w * \p x to the state)
  void accumul(const double & w, const State & x);

  /// ATLAS FieldSet interface
  /// For models that are not using ATLAS fieldsets for their own State data:
  /// - "toFieldSet" allocates the ATLAS fieldset based on the variables present in the State and
  ///    copies State data into the fieldset, including halo.
  /// - "fromFieldSet" copies fieldset data back into the State (interior points only).
  /// For models that are using ATLAS fieldsets for their own Incerment data, fields are shared from
  /// a fieldset to another. A working example is available with the QUENCH testbed of SABER.
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

  /// Serialize and deserialize (used in 4DEnVar, weak-constraint 4DVar and Block-Lanczos minimizer)
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  std::unique_ptr<State_> state_;
  size_t ID_;
  void print(std::ostream &) const override;

 protected:
  mutable std::unique_ptr<FieldSet3D> fset_;
};

// =============================================================================

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const Variables & vars,
                    const util::DateTime & time)
  : state_(), ID_(0)
{
  Log::trace() << "State<MODEL>::State starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(resol.geometry(), vars, time));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const eckit::Configuration & config)
  : state_(), ID_(config.getUnsigned("ID", 0))
{
  Log::trace() << "State<MODEL>::State read starting" << std::endl;
  util::Timer timer(classname(), "State");
  ASSERT(ID_ >= 0);
  state_.reset(new State_(resol.geometry(), config));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const State & other)
  : state_(), ID_(0)
{
  Log::trace() << "State<MODEL>::State interpolated starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(resol.geometry(), *other.state_));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Variables & vars, const State & other)
  : state_(), ID_(0)
{
  Log::trace() << "State<MODEL>::State variables starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(vars, *other.state_));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State variables done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const State & other)
  : state_(), ID_(other.ID())
{
  Log::trace() << "State<MODEL>::State starting copy" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(*other.state_));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::~State() {
  Log::trace() << "State<MODEL>::~State starting" << std::endl;
  util::Timer timer(classname(), "~State");
  state_.reset();
  fset_.reset();
  Log::trace() << "State<MODEL>::~State done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> & State<MODEL>::operator=(const State & rhs) {
  Log::trace() << "State<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  if (fset_) fset_->clear();
  *state_ = *rhs.state_;
  Log::trace() << "State<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "State<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  if (ID_ == 0) {
    ID_ = conf.getUnsigned("ID", 0);
    ASSERT(ID_ >= 0);
  }
  state_->read(conf);
  Log::trace() << "State<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "State<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  state_->write(conf);
  Log::trace() << "State<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double State<MODEL>::norm() const {
  Log::trace() << "State<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = state_->norm();
  Log::trace() << "State<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const Variables & State<MODEL>::variables() const {
  Log::trace() << "State<MODEL>::variables starting" << std::endl;
  util::Timer timer(classname(), "variables");
  return state_->variables();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
size_t State<MODEL>::serialSize() const {
  Log::trace() << "State<MODEL>::serialSize" << std::endl;
  util::Timer timer(classname(), "serialSize");
  return state_->serialSize();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::serialize(std::vector<double> & vect) const {
  Log::trace() << "State<MODEL>::serialize starting" << std::endl;
  util::Timer timer(classname(), "serialize");
  state_->serialize(vect);
  Log::trace() << "State<MODEL>::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::deserialize(const std::vector<double> & vect, size_t & current) {
  Log::trace() << "State<MODEL>::State deserialize starting" << std::endl;
  util::Timer timer(classname(), "deserialize");
  if (fset_) fset_->clear();
  state_->deserialize(vect, current);
  Log::trace() << "State<MODEL>::State deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::toFieldSet(atlas::FieldSet & fset) const {
  Log::trace() << "State<MODEL>::toFieldSet starting" << std::endl;
  util::Timer timer(classname(), "toFieldSet");
  ASSERT(fset.empty());
  state_->toFieldSet(fset);
  Log::trace() << "State<MODEL>::toFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------


template<typename MODEL>
void State<MODEL>::fromFieldSet(const atlas::FieldSet & fset) {
  Log::trace() << "State<MODEL>::fromFieldSet starting" << std::endl;
  util::Timer timer(classname(), "fromFieldSet");
  state_->fromFieldSet(fset);
  Log::trace() << "State<MODEL>::fromFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::print(std::ostream & os) const {
  Log::trace() << "State<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *state_;
  Log::trace() << "State<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::zero() {
  Log::trace() << "State<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  if (fset_) fset_->clear();
  state_->zero();
  Log::trace() << "State<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::accumul(const double & zz, const State & xx) {
  Log::trace() << "State<MODEL>::accumul starting" << std::endl;
  util::Timer timer(classname(), "accumul");
  if (fset_) fset_->clear();
  state_->accumul(zz, *xx.state_);
  Log::trace() << "State<MODEL>::accumul done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_STATE_H_
