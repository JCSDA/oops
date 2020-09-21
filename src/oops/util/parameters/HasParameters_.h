/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_HASPARAMETERS__H_
#define OOPS_UTIL_PARAMETERS_HASPARAMETERS__H_

#include <type_traits>

#include "oops/util/parameters/Parameters.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Detect if \c T defines \c Parameters_ to be the name of a class derived from
/// oops::Parameters using public inheritance.
///
/// If it does, HasParameters_ will be defined as std::true_type, otherwise as std::false_type.
///
/// \note Here's how this works. Note first that the primary template is followed by its
/// specialization. Suppose `HasParameters_ `is instantiated with a single parameter `T`. If `T`
/// doesn't define a type called `Parameters_`, substitution of `T` into the specialisation will
/// fail, so the primary template will be used. If `T` does define `Parameters_`, but not as the
/// name of a subclass of the `Parameters` class from `oops`, then `T::Parameters_*` won't be
/// convertible to `Parameters*`, so `std::is_convertible<...>::value` will be false. As a result,
/// `enable_if` won't define a `type` typedef, so again a substitution failure will occur and the
/// primary template will be used instead of the specialisation. Otherwise, if `T` defines
/// `Parameters_` as the name of a subclass of `oops::Parameters`, the substitution of `T` into
/// `std::enable_if<...>` will succeed and `std::enable_if<...>::type` will be defined as `void`.
/// This is identical to the value of the second parameter in the primary template. With both
/// parameters of the specialization matching that of the primary template, the compiler will
/// choose the specialization over the primary template.

// Primary template
template<typename T, typename = void>
struct HasParameters_ : std::false_type {};

// Specialization
template<typename T>
struct HasParameters_<T,
                      typename std::enable_if<std::is_convertible<typename T::Parameters_*,
                                                                  Parameters*>::value>::type>
    : std::true_type {};

// -----------------------------------------------------------------------------

template<typename T, typename FallbackType, bool THasParameters_>
struct TParameters_IfAvailableElseFallbackType;

template<typename T, typename FallbackType>
struct TParameters_IfAvailableElseFallbackType<T, FallbackType, false> {
  typedef FallbackType type;
};

template<typename T, typename FallbackType>
struct TParameters_IfAvailableElseFallbackType<T, FallbackType, true> {
  typedef typename T::Parameters_ type;
};

/// \brief Resolved to \c T::Parameters_ if \c T defines \c Parameters_ to be the name of a class
/// derived from oops::Parameters using public inheritance; otherwise resolved to \c FallbackType.
template<typename T, typename FallbackType>
using TParameters_IfAvailableElseFallbackType_t =
typename TParameters_IfAvailableElseFallbackType<T, FallbackType, HasParameters_<T>::value>::type;

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_HASPARAMETERS__H_
