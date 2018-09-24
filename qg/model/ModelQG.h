/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_MODELQG_H_
#define QG_MODEL_MODELQG_H_

#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "model/GeometryQG.h"
#include "model/QgFortran.h"
#include "model/QgTraits.h"
#include "oops/base/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class ModelBias;
  class FieldsQG;
  class StateQG;

// -----------------------------------------------------------------------------
/// QG model definition.
/*!
 *  QG nonlinear model definition and configuration parameters.
 */

class ModelQG: public oops::ModelBase<QgTraits>,
               private util::ObjectCounter<ModelQG> {
 public:
  static const std::string classname() {return "qg::ModelQG";}

  ModelQG(const GeometryQG &, const eckit::Configuration &);
  ~ModelQG();

/// Prepare model integration
  void initialize(StateQG &) const;

/// Model integration
  void step(StateQG &, const ModelBias &) const;
  int saveTrajectory(StateQG &, const ModelBias &) const;

/// Finish model integration
  void finalize(StateQG &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  util::Duration tstep_;
  const GeometryQG geom_;
  const oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_MODELQG_H_
