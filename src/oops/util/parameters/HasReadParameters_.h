/*
 * (C) Crown copyright 2021, Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_HASREADPARAMETERS__H_
#define OOPS_UTIL_PARAMETERS_HASREADPARAMETERS__H_

#include <type_traits>

#include "oops/util/parameters/Parameters.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Detect if \c T defines \c ReadParameters_ to be the name of a class derived from
/// oops::Parameters using public inheritance.
///
/// If it does, HasReadParameters_ will be defined as std::true_type, otherwise as std::false_type.
template<typename T, typename = void>
struct HasReadParameters_ : std::false_type {};

template<typename T>
struct HasReadParameters_<T,
                          typename std::enable_if<std::is_convertible<typename T::ReadParameters_*,
                                                                      Parameters*>::value>::type>
    : std::true_type {};

// -----------------------------------------------------------------------------

template<typename T, typename FallbackType, typename Enable = void>
struct TReadParameters_IfAvailableElseFallbackType{
  typedef FallbackType type;
};

template<typename T, typename FallbackType>
struct TReadParameters_IfAvailableElseFallbackType<
    T, FallbackType, typename std::enable_if<HasReadParameters_<T>::value>::type> {
  typedef typename T::ReadParameters_ type;
};

/// \brief Resolved to \c T::ReadParameters_ if \c T defines \c ReadParameters_ to be the name of
/// a class derived from oops::Parameters using public inheritance; otherwise resolved to \c
/// FallbackType.
template<typename T, typename FallbackType>
using TReadParameters_IfAvailableElseFallbackType_t =
typename TReadParameters_IfAvailableElseFallbackType<T, FallbackType>::type;

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_HASREADPARAMETERS__H_
