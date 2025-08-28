/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_TLMQG_H_
#define QG_MODEL_TLMQG_H_

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "model/ModelQG.h"
#include "model/QgFortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class IncrementQG;
  class ModelBias;
  class ModelBiasIncrement;
// -----------------------------------------------------------------------------
/// QG linear model definition.
/*!
 *  QG linear model definition and configuration parameters.
 */

class TlmQG: public util::Printable,
             private util::ObjectCounter<TlmQG> {
 public:
  static const std::string classname() {return "qg::TlmQG";}
  static std::vector<std::string> names() {return {"QgTLM"};}

  TlmQG(const GeometryQG &, const eckit::Configuration &);
  ~TlmQG();

  /// Prepare model integration
  void initializeTL(IncrementQG &) const;
  void initializeAD(IncrementQG &) const;

  /// Model integration
  void stepTL(IncrementQG &, const ModelBiasIncrement &) const;
  void stepAD(IncrementQG &, ModelBiasIncrement &) const;
  void setTrajectory(const StateQG &, StateQG &, const ModelBias &);

  /// Finish model integration
  void finalizeTL(IncrementQG &) const;
  void finalizeAD(IncrementQG &) const;

/// Other utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return steptraj_;}

 private:
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, int >::iterator trajIter;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;

// Data
  F90model keyConfig_;
  util::Duration tstep_;
  const GeometryQG resol_;
  std::map< util::DateTime, F90flds> traj_;
  const util::Duration steptraj_;
  const ModelQG lrmodel_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_TLMQG_H_
