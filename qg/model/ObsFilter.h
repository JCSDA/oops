/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_OBSFILTER_H_
#define QG_MODEL_OBSFILTER_H_

#include <memory>
#include <ostream>

#include "eckit/config/LocalConfiguration.h"

#include "model/QgTraits.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace qg {
  class GomQG;
  template <typename DATATYPE> class ObsDataQG;
  class ObsDiagsQG;
  class ObsSpaceQG;
  class ObsVecQG;

/// @brief ObsFilter is an optional GeoVaLs saver
class ObsFilter : public util::Printable {
 public:
  ObsFilter(const ObsSpaceQG &, const eckit::Configuration &,
            std::shared_ptr<ObsDataQG<int> >, std::shared_ptr<ObsDataQG<float> >,
            const int iteration = 0);

  void preProcess() {}
  void priorFilter(const GomQG &);
  void postFilter(const GomQG &,
                  const ObsVecQG &,
                  const ObsVecQG &,
                  const ObsDiagsQG &) {}

  oops::Variables requiredVars() const {return novars_;}
  oops::ObsVariables requiredHdiagnostics() const {return noobsvars_;}

 private:
  void print(std::ostream &) const override {}
  const oops::Variables novars_;
  const oops::ObsVariables noobsvars_;
  const eckit::LocalConfiguration config_;
  const bool saveGeoVaLs_;
};

}  // namespace qg

#endif  // QG_MODEL_OBSFILTER_H_
