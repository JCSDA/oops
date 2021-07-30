/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_OBSFILTERPARAMETERSBASE_H_
#define OOPS_GENERIC_OBSFILTERPARAMETERSBASE_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Base class of classes storing parameters controlling specific observation filters.
class ObsFilterParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(ObsFilterParametersBase, Parameters)
 public:
  /// \brief Observation filter type.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when observation filter parameters are deserialized into an
  /// ObsFilterParametersWrapper and used by FilterFactory to instantiate a filter whose type is
  /// determined at runtime), but not others (e.g. in tests written with a particular filter in
  /// mind). ObsFilterParametersWrapper will throw an exception if this parameter is not provided.
  OptionalParameter<std::string> filter{"filter", this};
};

}  // namespace oops

#endif  // OOPS_GENERIC_OBSFILTERPARAMETERSBASE_H_
