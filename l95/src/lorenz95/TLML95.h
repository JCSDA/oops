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
#include <vector>

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "lorenz95/ModelL95.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace lorenz95 {
  class FieldL95;
  class ModelTrajectory;
  class IncrementL95;
  class ModelBiasCorrection;

// -----------------------------------------------------------------------------
/// Lorenz 95 linear model definition.

class TLML95: public util::Printable,
              private util::ObjectCounter<TLML95> {
 public:
  static const std::string classname() {return "lorenz95::TLML95";}
  static std::vector<std::string> names() {return {"L95TLM"};}

  TLML95(const Resolution &, const eckit::Configuration &);
  ~TLML95();

/// Model trajectory computation
  void setTrajectory(const StateL95 &, StateL95 &, const ModelBias &);

/// Run TLM and its adjoint
  void initializeTL(IncrementL95 &) const;
  void stepTL(IncrementL95 &, const ModelBiasCorrection &) const;
  void finalizeTL(IncrementL95 &) const;

  void initializeAD(IncrementL95 &) const;
  void stepAD(IncrementL95 &, ModelBiasCorrection &) const;
  void finalizeAD(IncrementL95 &) const;

/// Other utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return steptraj_;}

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
  const util::Duration steptraj_;
  const double dt_;
  std::map< util::DateTime, ModelTrajectory * > traj_;
  const ModelL95 lrmodel_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_TLML95_H_
