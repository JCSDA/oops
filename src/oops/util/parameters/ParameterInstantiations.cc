/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cstdint>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// Explicit instantiations of the parameter class templates for the most commonly used value
// types, matching the extern template declarations at the end of Parameter.h,
// OptionalParameter.h and RequiredParameter.h.

template class Parameter<bool>;
template class Parameter<int>;
template class Parameter<size_t>;
template class Parameter<int64_t>;
template class Parameter<float>;
template class Parameter<double>;
template class Parameter<std::string>;
template class Parameter<util::DateTime>;
template class Parameter<util::Duration>;
template class Parameter<eckit::LocalConfiguration>;
template class Parameter<std::vector<std::string>>;
template class Parameter<std::vector<int>>;
template class Parameter<std::vector<float>>;
template class Parameter<std::vector<double>>;

template class OptionalParameter<bool>;
template class OptionalParameter<int>;
template class OptionalParameter<size_t>;
template class OptionalParameter<int64_t>;
template class OptionalParameter<float>;
template class OptionalParameter<double>;
template class OptionalParameter<std::string>;
template class OptionalParameter<util::DateTime>;
template class OptionalParameter<util::Duration>;
template class OptionalParameter<eckit::LocalConfiguration>;
template class OptionalParameter<std::vector<std::string>>;
template class OptionalParameter<std::vector<int>>;
template class OptionalParameter<std::vector<float>>;
template class OptionalParameter<std::vector<double>>;

template class RequiredParameter<bool>;
template class RequiredParameter<int>;
template class RequiredParameter<size_t>;
template class RequiredParameter<int64_t>;
template class RequiredParameter<float>;
template class RequiredParameter<double>;
template class RequiredParameter<std::string>;
template class RequiredParameter<util::DateTime>;
template class RequiredParameter<util::Duration>;
template class RequiredParameter<eckit::LocalConfiguration>;
template class RequiredParameter<std::vector<std::string>>;
template class RequiredParameter<std::vector<int>>;
template class RequiredParameter<std::vector<float>>;
template class RequiredParameter<std::vector<double>>;

}  // namespace oops
