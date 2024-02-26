/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_MODELBIAS_H_
#define LORENZ95_MODELBIAS_H_

#include <cmath>
#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class Resolution;
  class ModelBiasCorrection;

// -----------------------------------------------------------------------------

/// Model error for Lorenz 95 model.
/*!
 * This class is used to manipulate parameters of the model that
 * can be estimated in the assimilation. This includes model bias for
 * example but could be used for other parameters to be estimated.
 */

// -----------------------------------------------------------------------------

class ModelBias : public util::Printable,
                  private boost::noncopyable,
                  private util::ObjectCounter<ModelBias> {
 public:
  static const std::string classname() {return "lorenz95::ModelBias";}

  ModelBias(const Resolution &, const eckit::Configuration &);
  ModelBias(const Resolution &, const ModelBias &);
  ModelBias(const ModelBias &, const bool);
  ~ModelBias() {}

  ModelBias & operator=(const ModelBias &);
  ModelBias & operator+=(const ModelBiasCorrection &);

  const double & bias() const {return bias_;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return std::abs(bias_);}

 private:
  void print(std::ostream &) const;
  double bias_;
  bool active_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_MODELBIAS_H_
