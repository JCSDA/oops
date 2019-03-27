/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_LOCALIZATION_H_
#define OOPS_BASE_LOCALIZATION_H_

#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/assimilation/Increment4D.h"
#include "oops/base/StateEnsemble.h"
#include "oops/generic/LocalizationGeneric.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/LocalizationBase.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class Localization : public util::Printable,
                     private boost::noncopyable,
                     private util::ObjectCounter<Localization<MODEL>> {
  typedef LocalizationBase<MODEL>                 LocalizationBase_;
  typedef LocalizationGeneric<MODEL>              LocalizationGeneric_;
  typedef Increment<MODEL>                        Increment_;
  typedef Increment4D<MODEL>                      Increment4D_;
  typedef Geometry<MODEL>                         Geometry_;
  typedef boost::shared_ptr<StateEnsemble<MODEL>> EnsemblePtr_;

 public:
  static const std::string classname() {return "oops::Localization";}

  Localization(const Geometry_ &,
               const EnsemblePtr_,
               const eckit::Configuration &);
  virtual ~Localization();

  void multiply(Increment_ &) const;
  void multiply(Increment4D_ &) const;

 private:
  void print(std::ostream &) const;
  unsigned nsubwin_;
  boost::ptr_vector<LocalizationBase_> locBase_;
  boost::scoped_ptr<LocalizationGeneric_> locGen_;
};

// =============================================================================

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & resol,
                                  const EnsemblePtr_ ens,
                                  const eckit::Configuration & conf) : locBase_(), locGen_()
{
  Log::trace() << "Localization<MODEL>::Localization starting" << std::endl;
  util::Timer timer(classname(), "Localization");

  locGen_.reset(LocalizationGenericFactory<MODEL>::create(resol, ens, conf));
  if (!locGen_) {
    if (conf.has("localization_time")) {
      // Specific localization for each time-slot
      std::vector<eckit::LocalConfiguration> confTime;
      conf.get("covariance_time", confTime);
      nsubwin_ = confTime.size();
      for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
        locBase_.push_back(LocalizationFactory<MODEL>::create(resol, confTime[jsub]));
      }
    } else {
      // Only one localization
     nsubwin_ = 1;
     locBase_.push_back(LocalizationFactory<MODEL>::create(resol, conf));
    }
  }
  Log::trace() << "Localization<MODEL>::Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  Log::trace() << "Localization<MODEL>::~Localization starting" << std::endl;
  util::Timer timer(classname(), "~Localization");
  if (locGen_) {
    locGen_.reset();
  } else {
    locBase_.clear();
  }
  Log::trace() << "Localization<MODEL>::~Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  Log::trace() << "Localization<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  if (locGen_) {
    locGen_->doMultiply(dx);
  } else {
    locBase_[0].doMultiply(dx);
  }
  Log::trace() << "Localization<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment4D_ & dx) const {
  Log::trace() << "Localization<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  if (locGen_) {
    locGen_->doMultiply(dx);
  } else {
    if (nsubwin_ == 1) {
      // Only one localization
      locBase_[0].doMultiply(dx);
    } else {
      // Specific localization for each time-slot
      ASSERT(locBase_.size()-(dx.last()-dx.first()+1) == 0);
      for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
        int isub = jsub+dx.first();
        locBase_[jsub].doMultiply(dx[isub]);
      }
    }
  }
  Log::trace() << "Localization<MODEL>::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Localization<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  if (locGen_) {
    os << *locGen_;
  } else {
    for (unsigned jsub = 0; jsub < nsubwin_; ++jsub) {
      os << locBase_[jsub];
    }
  }
  Log::trace() << "Localization<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_LOCALIZATION_H_
