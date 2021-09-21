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

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/ModelQG.h"
#include "oops/qg/QgFortran.h"
#include "oops/qg/QgTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
// -----------------------------------------------------------------------------
/// QG linear model definition.
/*!
 *  QG linear model definition and configuration parameters.
 */

class TlmQG: public oops::interface::LinearModelBase<QgTraits>,
             private util::ObjectCounter<TlmQG> {
 public:
  static const std::string classname() {return "qg::TlmQG";}

  TlmQG(const GeometryQG &, const eckit::Configuration &);
  ~TlmQG();

  /// Prepare model integration
  void initializeTL(IncrementQG &) const override;
  void initializeAD(IncrementQG &) const override;

  /// Model integration
  void stepTL(IncrementQG &, const ModelBiasIncrement &) const override;
  void stepAD(IncrementQG &, ModelBiasIncrement &) const override;
  void setTrajectory(const StateQG &, StateQG &, const ModelBias &) override;

  /// Finish model integration
  void finalizeTL(IncrementQG &) const override;
  void finalizeAD(IncrementQG &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const GeometryQG & resolution() const {return resol_;}
  const oops::Variables & variables() const override {return linvars_;}

 private:
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, int >::iterator trajIter;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;

// Data
  F90model keyConfig_;
  util::Duration tstep_;
  const GeometryQG resol_;
  std::map< util::DateTime, F90flds> traj_;
  const ModelQG lrmodel_;
  oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_TLMQG_H_
