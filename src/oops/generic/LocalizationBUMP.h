/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_LOCALIZATIONBUMP_H_
#define OOPS_GENERIC_LOCALIZATIONBUMP_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/generic/LocalizationGeneric.h"
#include "oops/generic/OoBump.h"
#include "oops/generic/ParametersBUMP.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// BUMP localization matrix.

template<typename MODEL> class LocalizationBUMP : public LocalizationGeneric<MODEL> {
  typedef Geometry<MODEL>                         Geometry_;
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef OoBump<MODEL>                           OoBump_;
  typedef ParametersBUMP<MODEL>                   Parameters_;
  typedef IncrementEnsemble<MODEL>                Ensemble_;
  typedef boost::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  LocalizationBUMP(const Geometry_ &,
                   const EnsemblePtr_,
                   const eckit::Configuration &);
  ~LocalizationBUMP();

  void multiply(Increment_ &) const;
  void multiply(Increment4D_ &) const;

 private:
  void print(std::ostream &) const;

  std::unique_ptr<OoBump_> ooBump_;
  std::vector<util::DateTime> timeslots_;
};

// =============================================================================

template<typename MODEL>
LocalizationBUMP<MODEL>::LocalizationBUMP(const Geometry_ & resol,
                                          const EnsemblePtr_ ens,
                                          const eckit::Configuration & conf)
  : ooBump_()
{
// Setup variables
  const Variables vars(conf);

//  Setup timeslots
  const std::vector<std::string> timeslots_str(conf.getStringVector("timeslots"));
  std::vector<util::DateTime> timeslots;
  for (unsigned jsub = 0; jsub < timeslots_str.size(); ++jsub) {
    const util::DateTime date(timeslots_str[jsub]);
    timeslots.push_back(date);
  }
  Log::info() << "Number of ensemble time-slots:" << timeslots.size() << std::endl;

// Setup parameters
  Parameters_ param(resol, vars, timeslots, conf, ens);

// Transfer OoBump pointer
  ooBump_.reset(new OoBump_(param.getOoBump()));

  Log::trace() << "LocalizationBUMP:LocalizationBUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationBUMP<MODEL>::~LocalizationBUMP() {
  Log::trace() << "LocalizationBUMP:~LocalizationBUMP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationBUMP:multiply starting" << std::endl;
  ooBump_->multiplyNicas(dx);
  Log::trace() << "LocalizationBUMP:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::multiply(Increment4D_ & dx) const {
  Log::trace() << "LocalizationBUMP:multiply starting" << std::endl;
  ooBump_->multiplyNicas(dx);
  Log::trace() << "LocalizationBUMP:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::print(std::ostream & os) const {
  os << "LocalizationBUMP<MODEL>::print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONBUMP_H_
