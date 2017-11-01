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
#include <string>

#include "util/DateTime.h"

namespace eckit {
  class Configuration;
}

namespace util {
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
  virtual ~WeightingFct() {}

  virtual std::map< util::DateTime, double > setWeights(const util::DateTime &,
                                                        const util::DateTime &,
                                                        const util::Duration &) =0;
};

// -----------------------------------------------------------------------------

/// Factory
class WeightFactory {
 public:
  static WeightingFct * create(const eckit::Configuration &);
  virtual ~WeightFactory() { makers_->clear(); }

 protected:
  explicit WeightFactory(const std::string &);
 private:
  virtual WeightingFct * make(const eckit::Configuration &) =0;
  static std::map < std::string, WeightFactory * > * makers_;
};

template<class FCT>
class WeightMaker : public WeightFactory {
  virtual WeightingFct * make(const eckit::Configuration & config)
    {return new FCT(config);}
 public:
  explicit WeightMaker(const std::string & name) : WeightFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WEIGHTINGFCT_H_
