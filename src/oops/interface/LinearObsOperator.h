/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LINEAROBSOPERATOR_H_
#define OOPS_INTERFACE_LINEAROBSOPERATOR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename OBS>
class LinearObsOperator : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<LinearObsOperator<OBS> > {
  typedef typename OBS::LinearObsOperator  LinearObsOper_;
  typedef GeoVaLs<OBS>             GeoVaLs_;
  typedef ObsAuxControl<OBS>       ObsAuxControl_;
  typedef ObsAuxIncrement<OBS>     ObsAuxIncrement_;
  typedef ObsSpace<OBS>            ObsSpace_;
  typedef ObsVector<OBS>           ObsVector_;

 public:
  static const std::string classname() {return "oops::LinearObsOperator";}

  LinearObsOperator(const ObsSpace_ &, const eckit::Configuration &);
  ~LinearObsOperator();

/// Interfacing
  const LinearObsOper_ & linearobsoperator() const {return *oper_;}

/// Obs Operators
  void setTrajectory(const GeoVaLs_ &, const ObsAuxControl_ &);
  void simulateObsTL(const GeoVaLs_ &, ObsVector_ &, const ObsAuxIncrement_ &) const;
  void simulateObsAD(GeoVaLs_ &, const ObsVector_ &, ObsAuxIncrement_ &) const;

/// Other
  const Variables & requiredVars() const;  // Required inputs variables from LinearModel

 private:
  void print(std::ostream &) const;
  std::unique_ptr<LinearObsOper_> oper_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
LinearObsOperator<OBS>::LinearObsOperator(const ObsSpace_ & os,
                                            const eckit::Configuration & config): oper_() {
  Log::trace() << "LinearObsOperator<OBS>::LinearObsOperator starting" << std::endl;
  util::Timer timer(classname(), "LinearObsOperator");
  oper_.reset(new LinearObsOper_(os.obsspace(), config));
  Log::trace() << "LinearObsOperator<OBS>::LinearObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
LinearObsOperator<OBS>::~LinearObsOperator() {
  Log::trace() << "LinearObsOperator<OBS>::~LinearObsOperator starting" << std::endl;
  util::Timer timer(classname(), "~LinearObsOperator");
  oper_.reset();
  Log::trace() << "LinearObsOperator<OBS>::~LinearObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void LinearObsOperator<OBS>::setTrajectory(const GeoVaLs_ & gvals, const ObsAuxControl_ & aux) {
  Log::trace() << "LinearObsOperator<OBS>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  oper_->setTrajectory(gvals.geovals(), aux.obsauxcontrol());
  Log::trace() << "LinearObsOperator<OBS>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void LinearObsOperator<OBS>::simulateObsTL(const GeoVaLs_ & gvals, ObsVector_ & yy,
                                             const ObsAuxIncrement_ & aux) const {
  Log::trace() << "LinearObsOperator<OBS>::simulateObsTL starting" << std::endl;
  util::Timer timer(classname(), "simulateObsTL");
  oper_->simulateObsTL(gvals.geovals(), yy.obsvector(), aux.obsauxincrement());
  Log::trace() << "LinearObsOperator<OBS>::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void LinearObsOperator<OBS>::simulateObsAD(GeoVaLs_ & gvals, const ObsVector_ & yy,
                                             ObsAuxIncrement_ & aux) const {
  Log::trace() << "LinearObsOperator<OBS>::simulateObsAD starting" << std::endl;
  util::Timer timer(classname(), "simulateObsAD");
  oper_->simulateObsAD(gvals.geovals(), yy.obsvector(), aux.obsauxincrement());
  Log::trace() << "LinearObsOperator<OBS>::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
const Variables & LinearObsOperator<OBS>::requiredVars() const {
  Log::trace() << "LinearObsOperator<OBS>::requiredVars starting" << std::endl;
  util::Timer timer(classname(), "requiredVars");
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

template<typename OBS>
void LinearObsOperator<OBS>::print(std::ostream & os) const {
  Log::trace() << "LinearObsOperator<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *oper_;
  Log::trace() << "LinearObsOperator<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEAROBSOPERATOR_H_
