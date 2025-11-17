/*
 * (C) Copyright 2020-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/ObsVector.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Printable.h"

namespace oops {

template<typename MODEL, typename OBS>
class ObsLocalization : public util::Printable,
                        private boost::noncopyable {
  typedef typename MODEL::ObsLocalization     ObsLocalization_;
  typedef GeometryIterator<MODEL>    Iterator_;
  typedef ObsVector<OBS>             ObsVector_;
  typedef ObsSpace<OBS>              ObsSpace_;

 public:
  ObsLocalization(const eckit::Configuration &, const ObsSpace_ &);
  ~ObsLocalization();

  void computeLocalization(const Iterator_ &, ObsVector_ &) const;

 private:
  void print(std::ostream &) const override;

  std::string name_;
  std::unique_ptr<ObsLocalization_> obsloc_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ObsLocalization<MODEL, OBS>::ObsLocalization(const eckit::Configuration & config,
                                             const ObsSpace_ & os)
  : name_("oops::ObsLocalization::"+os.obsname()), obsloc_()
{
  Log::trace() << "ObsLocalization<OBS>::ObsLocalization start" << std::endl;
  util::Timer timer(name_, "ObsLocalization");
  obsloc_.reset(new ObsLocalization_(config, os.obsspace()));
  Log::trace() << "ObsLocalization<OBS>::ObsLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ObsLocalization<MODEL, OBS>::~ObsLocalization() {
  Log::trace() << "ObsLocalization<OBS>::~ObsLocalization start" << std::endl;
  util::Timer timer(name_, "~ObsLocalization");
  obsloc_.reset();
  Log::trace() << "ObsLocalization<OBS>::~ObsLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void ObsLocalization<MODEL, OBS>::computeLocalization(const Iterator_ & giter,
                                                      ObsVector_ & loc) const {
  Log::trace() << "ObsLocalization<OBS>::computeLocalization start" << std::endl;
  util::Timer timer(name_, "computeLocalization");
  obsloc_->computeLocalization(giter.geometryiter(), loc.obsvector());
  Log::trace() << "ObsLocalization<OBS>::computeLocalization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void ObsLocalization<MODEL, OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsLocalization<OBS>::print starting" << std::endl;
  util::Timer timer(name_, "print");
  os << *obsloc_;
  Log::trace() << "ObsLocalization<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
