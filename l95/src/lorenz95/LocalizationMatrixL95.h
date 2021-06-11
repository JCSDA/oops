/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_LOCALIZATIONMATRIXL95_H_
#define LORENZ95_LOCALIZATIONMATRIXL95_H_

#include <ostream>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "lorenz95/L95Traits.h"

// Forward declarations
namespace lorenz95 {
  class IncrementL95;

/// Localization matrix for Lorenz 95 model.

// -----------------------------------------------------------------------------
class LocalizationMatrixL95: public util::Printable,
                             private util::ObjectCounter<LocalizationMatrixL95> {
 public:
  static const std::string classname() {return "lorenz95::LocalizationMatrixL95";}

  LocalizationMatrixL95(const Resolution &, const eckit::Configuration &);
  void randomize(IncrementL95 &) const;
  void multiply(IncrementL95 &) const;

 private:
  void print(std::ostream &) const override;
  const unsigned int resol_;
  const double rscale_;
  std::vector<double> coefs_;
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_LOCALIZATIONMATRIXL95_H_
