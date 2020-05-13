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

#ifndef QG_MODEL_LOCALIZATIONMATRIXQG_H_
#define QG_MODEL_LOCALIZATIONMATRIXQG_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "model/GeometryQG.h"
#include "model/QgFortran.h"
#include "model/QgTraits.h"

// Forward declarations
namespace qg {
  class GeometryQG;
  class IncrementQG;

/// Localization matrix for QG model.

// -----------------------------------------------------------------------------
class LocalizationMatrixQG: public util::Printable,
                            private util::ObjectCounter<LocalizationMatrixQG> {
 public:
  static const std::string classname() {return "qg::LocalizationMatrixQG";}

  LocalizationMatrixQG(const GeometryQG &, const eckit::Configuration &);
  ~LocalizationMatrixQG();

  void multiply(IncrementQG &) const;

 private:
  void print(std::ostream &) const override;
  F90lclz keyLocal_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_LOCALIZATIONMATRIXQG_H_
