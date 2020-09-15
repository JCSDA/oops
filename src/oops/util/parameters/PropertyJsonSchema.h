/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_PROPERTYJSONSCHEMA_H_
#define OOPS_UTIL_PARAMETERS_PROPERTYJSONSCHEMA_H_

#include <map>
#include <string>

namespace oops {

/// \brief A representation of the JSON schema used to validate a property of a JSON node of type
/// 'object'.
///
/// Example: the following object represents a JSON schema constraining a property to be
/// integer-valued and positive:
/// \code
/// PropertyJsonSchema schema{{"type", "integer"}, {"minimum", "1"}};
/// \endcode
typedef std::map<std::string, std::string> PropertyJsonSchema;

/// \brief Return a string containing the JSON schema represented by \p schema.
std::string toString(const PropertyJsonSchema &schema);

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PROPERTYJSONSCHEMA_H_
