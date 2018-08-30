/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_ESTIMATEPARAMS_H_
#define OOPS_RUNS_ESTIMATEPARAMS_H_

#include <sstream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/generic/ParametersBUMP.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EstimateParams : public Application {
  typedef ParametersBUMP<MODEL>        Parameters_;

 public:
// -----------------------------------------------------------------------------
  static const std::string classname() {return "oops::EstimateParams";}
  EstimateParams() {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~EstimateParams() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    util::Timer timer(classname(), "write");

    // Setup parameters
    Parameters_ param(fullConfig);

    // Estimate parameters
    param.estimate();

    // Write parameters
    param.write();

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::EstimateParams<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ESTIMATEPARAMS_H_
