/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_TLMIDQG_H_
#define QG_MODEL_TLMIDQG_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/LinearModelBase.h"

#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "model/QgTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
// -----------------------------------------------------------------------------
/// QG linear identity model definition.
/*!
 *  QG linear identity model definition and configuration parameters.
 */

class TlmIdQG: public oops::LinearModelBase<QgTraits>,
              private util::ObjectCounter<TlmIdQG> {
 public:
  static const std::string classname() {return "qg::TlmIdQG";}

  TlmIdQG(const GeometryQG &, const eckit::Configuration &);
  ~TlmIdQG();

/// Model trajectory computation
  void setTrajectory(const StateQG &, StateQG &, const ModelBias &) override;

/// Run TLM and its adjoint
  void initializeTL(IncrementQG &) const override;
  void stepTL(IncrementQG &, const ModelBiasIncrement &) const override;
  void finalizeTL(IncrementQG &) const override;

  void initializeAD(IncrementQG &) const override;
  void stepAD(IncrementQG &, ModelBiasIncrement &) const override;
  void finalizeAD(IncrementQG &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const GeometryQG & resolution() const {return resol_;}
  const oops::Variables & variables() const override {return linvars_;}

 private:
  void print(std::ostream &) const override;

// Data
  int keyConfig_;
  util::Duration tstep_;
  const GeometryQG resol_;
  const oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_TLMIDQG_H_
