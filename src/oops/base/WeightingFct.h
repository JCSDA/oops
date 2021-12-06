/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_WEIGHTINGFCT_H_
#define OOPS_BASE_WEIGHTINGFCT_H_

#include <map>

namespace util {
  class DateTime;
  class Duration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Weighting Function
/*!
 * Abstract base class for weighting functions for various filters.
 */

class WeightingFct {
 public:
  virtual ~WeightingFct() = default;

  virtual std::map< util::DateTime, double > setWeights(const util::DateTime &,
                                                        const util::DateTime &,
                                                        const util::Duration &) = 0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WEIGHTINGFCT_H_
