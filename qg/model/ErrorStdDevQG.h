/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_ERRORSTDDEVQG_H_
#define QG_MODEL_ERRORSTDDEVQG_H_


#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/GeometryQG.h"
#include "oops/qg/QgFortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GeometryQG;
  class StateQG;
  class IncrementQG;

// -----------------------------------------------------------------------------
/// QG background error standard deviations

class ErrorStdDevQG: public util::Printable {
 public:
  static const std::string classname() {return "qg::ErrorStdDevQG";}

  ErrorStdDevQG(const StateQG &, const StateQG &,
                const GeometryQG &, const eckit::Configuration &);
  ~ErrorStdDevQG();

/// Perform linear transforms
  void multiply(const IncrementQG &, IncrementQG &) const;
  void multiplyInverse(const IncrementQG &, IncrementQG &) const;
  void multiplyAD(const IncrementQG &, IncrementQG &) const;
  void multiplyInverseAD(const IncrementQG &, IncrementQG &) const;

 private:
  void print(std::ostream &) const override;

  F90error_stddev keyFtnConfig_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_ERRORSTDDEVQG_H_
