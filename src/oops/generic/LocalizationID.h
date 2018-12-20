/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_LOCALIZATIONID_H_
#define OOPS_GENERIC_LOCALIZATIONID_H_

#include <sstream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/EnsemblesCollection.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/LocalizationBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Identity localization matrix for fast testing

template<typename MODEL>
class LocalizationID : public LocalizationBase<MODEL> {
  typedef typename MODEL::Geometry   Geometry_;
  typedef typename MODEL::Increment  Increment_;

 public:
  LocalizationID(const Geometry_ &, const eckit::Configuration &);
  ~LocalizationID();

  void multiply(Increment_ &) const;

 private:
  void print(std::ostream &) const;
  int keyLID_;
};

// =============================================================================

template<typename MODEL>
LocalizationID<MODEL>::LocalizationID(const Geometry_ & resol,
                                      const eckit::Configuration & conf)
  : keyLID_(0)
{
  Log::trace() << "LocalizationID:LocalizationID constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationID<MODEL>::~LocalizationID() {
  Log::trace() << "LocalizationID:~LocalizationID destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationID<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "LocalizationID:multiply starting" << std::endl;
  Log::trace() << "LocalizationID:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationID<MODEL>::print(std::ostream & os) const {
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONID_H_
