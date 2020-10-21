/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_MODELBIASCOVARIANCE_H_
#define LORENZ95_MODELBIASCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {
  class ModelBias;
  class ModelBiasCorrection;
  class Resolution;

// -----------------------------------------------------------------------------

class ModelBiasCovariance : public util::Printable,
                            private boost::noncopyable,
                            private util::ObjectCounter<ModelBiasCovariance> {
 public:
  static const std::string classname() {return "lorenz95::ModelBiasCovariance";}

/// Constructor, destructor
  ModelBiasCovariance(const eckit::Configuration &, const Resolution &);
  ~ModelBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ModelBias &, const Resolution &) {}
  void multiply(const ModelBiasCorrection &, ModelBiasCorrection &) const;
  void inverseMultiply(const ModelBiasCorrection &, ModelBiasCorrection &) const;
  void randomize(ModelBiasCorrection &) const;

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream &) const;
  const eckit::LocalConfiguration conf_;
  double variance_;
  bool active_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_MODELBIASCOVARIANCE_H_
