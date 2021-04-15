/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_OBSLOCQG_H_
#define QG_MODEL_OBSLOCQG_H_

#include <ostream>

#include "atlas/field.h"

#include "oops/util/Printable.h"

#include "oops/qg/ObsDataQG.h"

namespace eckit {
  class Configuration;
}

namespace qg {
  class GeometryQGIterator;
  class ObsSpaceQG;
  class ObsVecQG;

/// \brief Observation-space localization for QG model (Heaviside function
/// with prescribed lengthscale).
class ObsLocQG : public util::Printable {
 public:
  ObsLocQG(const eckit::Configuration &, const ObsSpaceQG &);

  void computeLocalization(const GeometryQGIterator &, ObsDataQG<int> &, ObsVecQG &) const;

 private:
  void print(std::ostream &) const override;
  const double lengthscale_;
  const ObsSpaceQG & obsdb_;
};

}  // namespace qg

#endif  // QG_MODEL_OBSLOCQG_H_
