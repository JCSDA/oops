/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_STRINGFUNCTIONS_H_
#define OOPS_UTIL_STRINGFUNCTIONS_H_

#include <string>

#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"

namespace util {
  namespace stringfunctions {

    /// \brief Replace a placeholder with ensemble member index.
    ///
    /// If \p conf contains an integer-valued key "member", the "%{member}%" placeholder in \p
    /// filename is replaced by "XXX", where XXX is the value of \p member padded with zeros
    /// from the left to a total of \p ndigits digits.
    void swapNameMember(const eckit::Configuration &member, std::string &filename, int ndigits = 3);

    /// \brief Replace a placeholder with ensemble member index.
    ///
    /// \brief If \p member is not none, the "%{member}%" placeholder in \p filename is replaced
    /// by "XXX", where XXX is the value of \p member padded with zeros from the left
    /// to a total of \p ndigits digits.
    void swapNameMember(const boost::optional<int> &member, std::string &filename, int ndigits = 3);

    /// \brief Convert sequence elements to strings and join them using a delimiter.
    ///
    /// This is a generalization of eckit::StringTools::join. It calls \p toString on
    /// each element of the range [begin, end) and joins the resulting strings, inserting the
    /// delimiter \p delimiter between each pair of adjacent strings.
    ///
    /// \param toString
    ///   An unary function object taking a dereferenced \c Iterator and returning a \c std::string.
    template<typename Iterator, typename Converter>
    std::string join(const std::string &delimiter, Iterator begin, Iterator end,
                     Converter toString) {
      if (begin == end)
        return "";
      std::string r(toString(*begin));
      for (Iterator it = ++begin; it != end; ++it) {
        r += delimiter;
        r += toString(*it);
      }
      return r;
    }

  }  // namespace stringfunctions
}  // namespace util

#endif  // OOPS_UTIL_STRINGFUNCTIONS_H_
