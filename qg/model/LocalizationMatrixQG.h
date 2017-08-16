/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_LOCALIZATIONMATRIXQG_H_
#define QG_MODEL_LOCALIZATIONMATRIXQG_H_

#include <ostream>
#include <string>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/LocalizationBase.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"

#include "model/QgFortran.h"
#include "model/QgTraits.h"

// Forward declarations
namespace qg {
  class IncrementQG;
  class StateQG;

/// Localization matrix for QG model.

// -----------------------------------------------------------------------------
class LocalizationMatrixQG: public oops::LocalizationBase<QgTraits>,
                            private util::ObjectCounter<LocalizationMatrixQG> {
 public:
  static const std::string classname() {return "qg::LocalizationMatrixQG";}

  LocalizationMatrixQG(const StateQG &, const eckit::Configuration &);
  ~LocalizationMatrixQG();
  void multiply(IncrementQG &) const;

 private:
  void print(std::ostream &) const;
  F90lclz keyFtnConfig_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_LOCALIZATIONMATRIXQG_H_
