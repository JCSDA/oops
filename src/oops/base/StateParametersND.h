/*
 * (C) Crown copyright 2021, Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATEPARAMETERSND_H_
#define OOPS_BASE_STATEPARAMETERSND_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/variant.hpp>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Parameters of a set of 3D states valid at different times, each used by a different member of
/// the MPI communicator in time.
template <typename MODEL>
class StateParameters4D : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateParameters4D, Parameters)

  typedef typename MODEL::State State_;

 public:
  /// Set to State_::Parameters_ if State_ provides a type called Parameters_ and to
  /// GenericParameters (a thin wrapper of an eckit::LocalConfiguration object) if not.
  typedef TParameters_IfAvailableElseFallbackType_t<State_, GenericParameters> StateParameters3D_;

  /// \brief Return the parameters of the 3D state to be used by the calling process.
  ///
  /// \param timeComm The MPI communicator in time.
  const StateParameters3D_ &at(const eckit::mpi::Comm &timeComm) const;

  /// Parameters of a set of 3D states valid at different times.
  RequiredParameter<std::vector<StateParameters3D_>> states{"states", this};
};

template <typename MODEL>
const typename StateParameters4D<MODEL>::StateParameters3D_ &StateParameters4D<MODEL>::at(
    const eckit::mpi::Comm &timeComm) const {
  ASSERT(states.value().size() == timeComm.size());
  return states.value()[timeComm.rank()];
}


/// Parameters of either a single 3D model state or a set of 3D states valid at different times,
/// each used by a different member of the MPI communicator in time.
template <typename MODEL>
class StateParametersND : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateParametersND, Parameters)

  typedef typename MODEL::State State_;

 public:
  /// Set to State_::Parameters_ if State_ provides a type called Parameters_ and to
  /// GenericParameters (a thin wrapper of an eckit::LocalConfiguration object) if not.
  typedef TParameters_IfAvailableElseFallbackType_t<State_, GenericParameters> StateParameters3D_;
  typedef StateParameters4D<MODEL> StateParameters4D_;
  typedef boost::variant<StateParameters3D_, StateParameters4D_> StateParametersVariant_;

  /// \brief The stored value.
  const StateParametersVariant_ &value() const { return value_; }

  /// \brief The stored value.
  operator const StateParametersVariant_ &() const { return value_; }

  /// \brief Parameters of the 3D state to be used by the calling process.
  ///
  /// \param timeComm The MPI communicator in time.
  const StateParameters3D_ &at(const eckit::mpi::Comm &comm) const;

  // Import all overloads of deserialize() from the base class. We will override the virtual one.
  using ParameterBase::deserialize;

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

 private:
  StateParametersVariant_ value_;
};

template <typename MODEL>
const typename StateParametersND<MODEL>::StateParameters3D_ &
StateParametersND<MODEL>::at(const eckit::mpi::Comm &comm) const {
  if (const StateParameters3D_* parameters3D = boost::get<StateParameters3D_>(&value_)) {
    return *parameters3D;
  } else if (const StateParameters4D_* parameters4D = boost::get<StateParameters4D_>(&value_)) {
    return parameters4D->at(comm);
  } else {
    throw eckit::BadValue("StateParametersND::at() called before deserialize()", Here());
  }
}

template <typename MODEL>
void StateParametersND<MODEL>::deserialize(util::CompositePath &path,
                                           const eckit::Configuration &config) {
  Parameters::deserialize(path, config);
  if (config.has("states")) {
    value_ = StateParameters4D_();
    boost::get<StateParameters4D_>(value_).deserialize(path, config);
  } else {
    value_ = StateParameters3D_();
    boost::get<StateParameters3D_>(value_).deserialize(path, config);
  }
}

template <typename MODEL>
void StateParametersND<MODEL>::serialize(eckit::LocalConfiguration &config) const {
  Parameters::serialize(config);
  if (const StateParameters3D_* parameters3D = boost::get<StateParameters3D_>(&value_)) {
    parameters3D->serialize(config);
  } else if (const StateParameters4D_* parameters4D = boost::get<StateParameters4D_>(&value_)) {
    parameters4D->serialize(config);
  }
}

template <typename MODEL>
ObjectJsonSchema StateParametersND<MODEL>::jsonSchema() const {
  StateParameters3D_ params3D;
  StateParameters4D_ params4D;
  ObjectJsonSchema schema3D = params3D.jsonSchema();
  ObjectJsonSchema schema4D = params4D.jsonSchema();

  std::vector<ConditionalObjectJsonSchema> allOf(1);
  // If the object has a 'states' property...
  allOf[0].if_   = ObjectJsonSchema({}, {"states"}, true);
  // then expect it to represent a set of 3D states to be used by different members of the MPI
  // communicator in time...
  allOf[0].then  = schema4D;
  // otherwise expect it to represent a single 3D state
  allOf[0].else_ = schema3D;

  ObjectJsonSchema mySchema(std::move(allOf));

  mySchema.combineWith(Parameters::jsonSchema());
  return mySchema;
}

}  // namespace oops

#endif  // OOPS_BASE_STATEPARAMETERSND_H_
