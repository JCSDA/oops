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

template <typename MODEL>
class LinearObsOperator : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<LinearObsOperator<MODEL> > {
  typedef typename MODEL::LinearObsOperator  LinearObsOper_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncrement_;
  typedef ObsSpace<MODEL>            ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

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

template <typename MODEL>
LinearObsOperator<MODEL>::LinearObsOperator(const ObsSpace_ & os,
                                            const eckit::Configuration & config): oper_() {
  Log::trace() << "LinearObsOperator<MODEL>::LinearObsOperator starting" << std::endl;
  util::Timer timer(classname(), "LinearObsOperator");
  oper_.reset(new LinearObsOper_(os.obsspace(), config));
  Log::trace() << "LinearObsOperator<MODEL>::LinearObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperator<MODEL>::~LinearObsOperator() {
  Log::trace() << "LinearObsOperator<MODEL>::~LinearObsOperator starting" << std::endl;
  util::Timer timer(classname(), "~LinearObsOperator");
  oper_.reset();
  Log::trace() << "LinearObsOperator<MODEL>::~LinearObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LinearObsOperator<MODEL>::setTrajectory(const GeoVaLs_ & gvals, const ObsAuxControl_ & aux) {
  Log::trace() << "LinearObsOperator<MODEL>::setTrajectory starting" << std::endl;
  util::Timer timer(classname(), "setTrajectory");
  oper_->setTrajectory(gvals.geovals(), aux.obsauxcontrol());
  Log::trace() << "LinearObsOperator<MODEL>::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LinearObsOperator<MODEL>::simulateObsTL(const GeoVaLs_ & gvals, ObsVector_ & yy,
                                             const ObsAuxIncrement_ & aux) const {
  Log::trace() << "LinearObsOperator<MODEL>::simulateObsTL starting" << std::endl;
  util::Timer timer(classname(), "simulateObsTL");
  oper_->simulateObsTL(gvals.geovals(), yy.obsvector(), aux.obsauxincrement());
  Log::trace() << "LinearObsOperator<MODEL>::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LinearObsOperator<MODEL>::simulateObsAD(GeoVaLs_ & gvals, const ObsVector_ & yy,
                                             ObsAuxIncrement_ & aux) const {
  Log::trace() << "LinearObsOperator<MODEL>::simulateObsAD starting" << std::endl;
  util::Timer timer(classname(), "simulateObsAD");
  oper_->simulateObsAD(gvals.geovals(), yy.obsvector(), aux.obsauxincrement());
  Log::trace() << "LinearObsOperator<MODEL>::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
const Variables & LinearObsOperator<MODEL>::requiredVars() const {
  Log::trace() << "LinearObsOperator<MODEL>::requiredVars starting" << std::endl;
  util::Timer timer(classname(), "requiredVars");
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearObsOperator<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LinearObsOperator<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *oper_;
  Log::trace() << "LinearObsOperator<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEAROBSOPERATOR_H_
