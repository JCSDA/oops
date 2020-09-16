/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_HASPARAMETERS__H_
#define OOPS_UTIL_PARAMETERS_HASPARAMETERS__H_

#include <type_traits>

#include <boost/type_traits/make_void.hpp>  // for void_t. Can be removed in C++17

namespace oops {

/// \brief Detect if \c T provides a type called \c Parameters_ or not.
///
/// If it does, HasParameters_ will be defined as std::true_type, otherwise as std::false_type.
template<typename T, typename = boost::void_t<>>
struct HasParameters_ : std::false_type {};

template<typename T>
struct HasParameters_<T, boost::void_t<typename T::Parameters_>> : std::true_type {};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_HASPARAMETERS__H_
