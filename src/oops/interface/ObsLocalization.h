/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_OBSLOCALIZATION_H_
#define OOPS_INTERFACE_OBSLOCALIZATION_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS, typename LOC>
class ObsLocalization : public ObsLocalizationBase<MODEL, OBS> {
  typedef GeometryIterator<MODEL>  GeometryIterator_;
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsDataVector<OBS, int>  ObsDataVector_;
  typedef ObsVector<OBS>           ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsLocalization";}

  ObsLocalization(const eckit::Configuration &, const ObsSpace_ &);
  ~ObsLocalization();

  void computeLocalization(const GeometryIterator_ &,
                           ObsDataVector_ &, ObsVector_ &) const override;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<LOC> obsloc_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS, typename LOC>
ObsLocalization<MODEL, OBS, LOC>::ObsLocalization(const eckit::Configuration & conf,
                                                  const ObsSpace_ & obspace)
  : obsloc_()
{
  Log::trace() << "ObsLocalization<MODEL, OBS, LOC>::ObsLocalization starting" << std::endl;
  util::Timer timer(classname(), "ObsLocalization");
  obsloc_.reset(new LOC(conf, obspace.obsspace()));
  Log::trace() << "ObsLocalization<MODEL, OBS, LOC>::ObsLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS, typename LOC>
ObsLocalization<MODEL, OBS, LOC>::~ObsLocalization() {
  Log::trace() << "ObsLocalization<MODEL, OBS, LOC>::~ObsLocalization starting" << std::endl;
  util::Timer timer(classname(), "~ObsLocalization");
  obsloc_.reset();
  Log::trace() << "ObsLocalization<MODEL, OBS, LOC>::~ObsLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS, typename LOC>
void ObsLocalization<MODEL, OBS, LOC>::computeLocalization(const GeometryIterator_ & p,
                                       ObsDataVector_ & local, ObsVector_ & obsvector) const {
  Log::trace() << "ObsLocalization<MODEL, OBS, LOC>:: computeLocalization starting" << std::endl;
  util::Timer timer(classname(), "computeLocalization");
  obsloc_->computeLocalization(p.geometryiter(), local.obsdatavector(), obsvector.obsvector());
  Log::trace() << "ObsLocalization<MODEL, OBS, LOC>:: computeLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS, typename LOC>
void ObsLocalization<MODEL, OBS, LOC>::print(std::ostream & os) const {
  os << *obsloc_;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSLOCALIZATION_H_
