/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_LOCALIZATIONID_H_
#define OOPS_GENERIC_LOCALIZATIONID_H_

#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/generic/LocalizationGeneric.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// Identity localization matrix for fast testing

template<typename MODEL>
class LocalizationID : public LocalizationGeneric<MODEL> {
  typedef Geometry<MODEL>                         Geometry_;
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef boost::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  LocalizationID(const Geometry_ &,
                 const EnsemblePtr_,
                 const eckit::Configuration &);
  ~LocalizationID();

  void multiply(Increment_ &) const;
  void multiply(Increment4D_ &) const;

 private:
  void print(std::ostream &) const;
  int cross_timeslot_;
};

// =============================================================================

template<typename MODEL>
LocalizationID<MODEL>::LocalizationID(const Geometry_ & resol,
                                      const EnsemblePtr_ ens,
                                      const eckit::Configuration & conf)
  : cross_timeslot_(0)
{
  if (conf.has("cross_timeslot")) {
    cross_timeslot_ = conf.getInt("cross_timeslot");
  }
  ASSERT(cross_timeslot_ == 0 || cross_timeslot_ == 1);
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
void LocalizationID<MODEL>::multiply(Increment4D_ & dx) const {
  Log::trace() << "LocalizationID:multiply starting" << std::endl;
  if (cross_timeslot_ == 1) {
    // Sum over timeslots
    Increment_ dxtmp(dx[dx.first()]);
    for (int isub = dx.first()+1; isub <= dx.last(); ++isub) {
       dxtmp.axpy(1.0, dx[isub], false);
    }

    // Copy result to all timeslots
    for (int isub = dx.first(); isub <= dx.last(); ++isub) {
       dx[isub].zero();
       dx[isub].axpy(1.0, dxtmp, false);
    }
  }
  Log::trace() << "LocalizationID:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationID<MODEL>::print(std::ostream & os) const {
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_LOCALIZATIONID_H_
