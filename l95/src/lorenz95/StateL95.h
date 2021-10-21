/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_STATEL95_H_
#define LORENZ95_STATEL95_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "lorenz95/FieldL95.h"
#include "lorenz95/Resolution.h"

#include "oops/base/Variables.h"
#include "oops/base/WriteParametersBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class IncrementL95;

// -----------------------------------------------------------------------------

/// \brief Parameters used to initialize a Lorenz95 model's state.
class StateL95Parameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StateL95Parameters, Parameters)

 public:
  /// \brief Validity date.
  oops::RequiredParameter<util::DateTime> date{"date", this};
  /// \brief File to load the state from. Either this option or `analytic_init` must be set.
  oops::OptionalParameter<std::string> filename{"filename", this};
  /// \brief Options used to generate the state on the fly. Either this option or `filename`
  /// must be set.
  oops::OptionalParameter<Field95GenerateParameters> analyticInit{"analytic init", this};
  /// \brief Ensemble member index.
  oops::OptionalParameter<int> member{"member", this};
};

// -----------------------------------------------------------------------------

/// \brief Parameters controlling the action of writing a Lorenz95 model's state to a file.
class StateL95WriteParameters : public oops::WriteParametersBase {
  OOPS_CONCRETE_PARAMETERS(StateL95WriteParameters, WriteParametersBase)

 public:
  oops::RequiredParameter<std::string> datadir{"datadir", this};
  oops::RequiredParameter<std::string> exp{"exp", this};
  oops::RequiredParameter<std::string> type{"type", this};
};

/// L95 model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class StateL95 : public util::Printable,
                 private util::ObjectCounter<StateL95> {
 public:
  typedef StateL95Parameters Parameters_;
  typedef StateL95WriteParameters WriteParameters_;

  static const std::string classname() {return "lorenz95::StateL95";}

/// Constructor, destructor
  StateL95(const Resolution &, const oops::Variables &, const util::DateTime &);
  StateL95(const Resolution &, const Parameters_ &);
  StateL95(const Resolution &, const StateL95 &);
  StateL95(const StateL95 &);
  virtual ~StateL95();
  StateL95 & operator=(const StateL95 &);

/// Interactions with increments
  StateL95 & operator+=(const IncrementL95 &);

// Utilities
  const FieldL95 & getField() const {return fld_;}
  FieldL95 & getField() {return fld_;}
  std::shared_ptr<const Resolution> geometry() const {
    std::shared_ptr<const Resolution> geom(new Resolution(fld_.resol()));
    return geom;
  }

  void read(const Parameters_ &);
  void write(const WriteParameters_ &) const;
  double norm () const {return fld_.rms();}
  const util::DateTime & validTime() const {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}
  util::DateTime & validTime() {return time_;}
  const oops::Variables & variables() const {return vars_;}

// For accumulator
  void zero();
  void accumul(const double &, const StateL95 &);

/// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

 private:
  void print(std::ostream &) const;
  FieldL95 fld_;
  util::DateTime time_;
  oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_STATEL95_H_
