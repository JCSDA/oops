/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSERVATIONTLAD_H_
#define LORENZ95_OBSERVATIONTLAD_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "lorenz95/ObservationL95.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace lorenz95 {
  class GomL95;
  class ObservationL95;
  class ObsBias;
  class ObsBiasCorrection;
  class ObsTable;
  class ObsVec1D;
  class NoVariables;

/// Observation for Lorenz 95 model.
/*!
 *  ObservationTLAD defines the linearized ObsOperator for Lorenz 95 model.
 */

// -----------------------------------------------------------------------------

class ObservationTLAD : public util::Printable,
                        private boost::noncopyable,
                        private util::ObjectCounter<ObservationTLAD> {
 public:
  static const std::string classname() {return "lorenz95::ObservationTLAD";}

  static ObservationTLAD * create(const ObservationL95 & hop)
    {return new ObservationTLAD(hop.table());}

  ~ObservationTLAD();

// Obs Operators
  void setTrajectory(const GomL95 &, const ObsBias &);
  void obsEquivTL(const GomL95 &, ObsVec1D &, const ObsBiasCorrection &) const;
  void obsEquivAD(GomL95 &, const ObsVec1D &, ObsBiasCorrection &) const;

// Other
  boost::shared_ptr<const NoVariables> variables() const {return inputs_;}

 private:
  void print(std::ostream &) const;
  explicit ObservationTLAD(const ObsTable &);
  boost::shared_ptr<const NoVariables> inputs_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSERVATIONTLAD_H_
