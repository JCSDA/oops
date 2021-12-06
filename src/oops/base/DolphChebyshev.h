/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_DOLPHCHEBYSHEV_H_
#define OOPS_BASE_DOLPHCHEBYSHEV_H_

#include <map>

#include "oops/base/WeightingFct.h"
#include "oops/util/Duration.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace oops {

// -----------------------------------------------------------------------------

class DolphChebyshev : public WeightingFct {
 public:
  explicit DolphChebyshev(const eckit::Configuration &);
  ~DolphChebyshev() {}

  std::map< util::DateTime, double > setWeights(const util::DateTime &,
                                                const util::DateTime &,
                                                const util::Duration &);

 private:
  util::Duration tau_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_DOLPHCHEBYSHEV_H_
