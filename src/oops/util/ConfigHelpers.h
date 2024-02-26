/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace eckit {
class LocalConfiguration;
}

// ConfigHelpers implement common OOPS/JEDI tasks with eckit Configurations; ConfigHelpers may
// call ConfigFunctions but not vice versa.

namespace util {

/// \brief Sets a member number for an ensemble, including performing the member template
///        pattern substitution if needed.
void setMember(eckit::LocalConfiguration & conf, int);

}  // namespace util
