/*
 * (C) Crown copyright 2021, Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <type_traits>

#include "oops/util/parameters/Parameters.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Detect if \c T defines \c DiracParameters_ to be the name of a class derived from
/// oops::Parameters using public inheritance.
///
/// If it does, HasDiracParameters_ will be defined as std::true_type, otherwise as std::false_type.
template<typename T, typename = void>
struct HasDiracParameters_ : std::false_type {};

template<typename T>
struct HasDiracParameters_<T,
                          typename std::enable_if<std::is_convertible<typename T::DiracParameters_*,
                                                                      Parameters*>::value>::type>
    : std::true_type {};

// -----------------------------------------------------------------------------

template<typename T, typename FallbackType, typename Enable = void>
struct TDiracParameters_IfAvailableElseFallbackType{
  typedef FallbackType type;
};

template<typename T, typename FallbackType>
struct TDiracParameters_IfAvailableElseFallbackType<
    T, FallbackType, typename std::enable_if<HasDiracParameters_<T>::value>::type> {
  typedef typename T::DiracParameters_ type;
};

/// \brief Resolved to \c T::DiracParameters_ if \c T defines \c DiracParameters_ to be the name of
/// a class derived from oops::Parameters using public inheritance; otherwise resolved to \c
/// FallbackType.
template<typename T, typename FallbackType>
using TDiracParameters_IfAvailableElseFallbackType_t =
typename TDiracParameters_IfAvailableElseFallbackType<T, FallbackType>::type;

}  // namespace oops

