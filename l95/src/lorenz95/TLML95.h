/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_TLML95_H_
#define LORENZ95_TLML95_H_

#include <map>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/interface/LinearModelBase.h"

#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "lorenz95/L95Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace lorenz95 {
  class FieldL95;

// -----------------------------------------------------------------------------
/// Lorenz 95 linear model definition.

class TLML95: public oops::LinearModelBase<L95Traits>,
              private util::ObjectCounter<TLML95> {
 public:
  static const std::string classname() {return "lorenz95::TLML95";}

  TLML95(const Resolution &, const eckit::Configuration &);
  ~TLML95();

/// Model trajectory computation
  void setTrajectory(const StateL95 &, StateL95 &, const ModelBias &) override;

/// Run TLM and its adjoint
  void initializeTL(IncrementL95 &) const override;
  void stepTL(IncrementL95 &, const ModelBiasCorrection &) const override;
  void finalizeTL(IncrementL95 &) const override;

  void initializeAD(IncrementL95 &) const override;
  void stepAD(IncrementL95 &, ModelBiasCorrection &) const override;
  void finalizeAD(IncrementL95 &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const Resolution & resolution() const {return resol_;}

 private:
  const ModelTrajectory * getTrajectory(const util::DateTime &) const;
  void tendenciesTL(const FieldL95 &, const double &, const FieldL95 &, FieldL95 &) const;
  void tendenciesAD(FieldL95 &, double &, const FieldL95 &, const FieldL95 &) const;
  void print(std::ostream &) const override;

  typedef std::map< util::DateTime, ModelTrajectory * >::iterator trajIter;
  typedef std::map< util::DateTime, ModelTrajectory * >::const_iterator trajICst;

// Data
  const Resolution resol_;
  const util::Duration tstep_;
  const double dt_;
  std::map< util::DateTime, ModelTrajectory * > traj_;
  const ModelL95 lrmodel_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_TLML95_H_
