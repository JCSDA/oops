/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSLOCALIZATION_H_
#define OOPS_INTERFACE_OBSLOCALIZATION_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename OBS, typename LOC>
class ObsLocalization : public ObsLocalizationBase<OBS> {
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsVector<OBS>           ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsLocalization";}

  ObsLocalization(const eckit::Configuration &, const ObsSpace_ &);
  ~ObsLocalization();

  void multiply(ObsVector_ &) const override;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<LOC> obsloc_;
};

// -----------------------------------------------------------------------------

template <typename OBS, typename LOC>
ObsLocalization<OBS, LOC>::ObsLocalization(const eckit::Configuration & conf,
                                             const ObsSpace_ & os)
  : obsloc_()
{
  Log::trace() << "ObsLocalization<OBS, LOC>::ObsLocalization Configuration starting" <<
                  std::endl;
  util::Timer timer(classname(), "ObsLocalization");
  obsloc_.reset(new LOC(conf, os.obsspace()));
  Log::trace() << "ObsLocalization<OBS, LOC>::ObsLocalization Configuration done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename LOC>
ObsLocalization<OBS, LOC>::~ObsLocalization() {
  Log::trace() << "ObsLocalization<OBS, LOC>::~ObsLocalization starting" << std::endl;
  util::Timer timer(classname(), "~ObsLocalization");
  obsloc_.reset();
  Log::trace() << "ObsLocalization<OBS, LOC>::~ObsLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename LOC>
void ObsLocalization<OBS, LOC>::multiply(ObsVector_ & dy) const {
  Log::trace() << "ObsLocalization<OBS, LOC>:: preProcess starting" << std::endl;
  util::Timer timer(classname(), "preProcess");
  obsloc_->multiply(dy.obsvector());
  Log::trace() << "ObsLocalization<OBS, LOC>:: preProcess done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS, typename LOC>
void ObsLocalization<OBS, LOC>::print(std::ostream & os) const {
  os << "ObsLocalization " << *obsloc_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSLOCALIZATION_H_
