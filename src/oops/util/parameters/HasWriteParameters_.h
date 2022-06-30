/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_HASWRITEPARAMETERS__H_
#define OOPS_UTIL_PARAMETERS_HASWRITEPARAMETERS__H_

#include <type_traits>

#include "oops/base/WriteParametersBase.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Detect if \c T defines \c WriteParameters_ to be the name of a class derived from
/// oops::WriteParametersBase using public inheritance.
///
/// If it does, HasWriteParameters_ will be defined as std::true_type, otherwise as std::false_type.
template<typename T, typename = void>
struct HasWriteParameters_ : std::false_type {};

template<typename T>
struct HasWriteParameters_<T,
                          typename std::enable_if<std::is_convertible<typename T::WriteParameters_*,
                                                                WriteParametersBase*>::value>::type>
    : std::true_type {};

// -----------------------------------------------------------------------------

template<typename T, typename FallbackType, typename Enable = void>
struct TWriteParameters_IfAvailableElseFallbackType{
  typedef FallbackType type;
};

template<typename T, typename FallbackType>
struct TWriteParameters_IfAvailableElseFallbackType<
    T, FallbackType, typename std::enable_if<HasWriteParameters_<T>::value>::type> {
  typedef typename T::WriteParameters_ type;
};

/// \brief Resolved to \c T::WriteParameters_ if \c T defines \c WriteParameters_ to be the name of
/// a class derived from oops::WriteParametersBase using public inheritance; otherwise resolved to
/// \c FallbackType.
template<typename T, typename FallbackType>
using TWriteParameters_IfAvailableElseFallbackType_t =
typename TWriteParameters_IfAvailableElseFallbackType<T, FallbackType>::type;

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_HASWRITEPARAMETERS__H_
