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

#include <string>

#include <boost/shared_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/Variables.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class LinearObsOperator : public util::Printable,
                          private util::ObjectCounter<LinearObsOperator<MODEL> > {
  typedef typename MODEL::LinearObsOperator     LinearObsOperator_;
  typedef ModelAtLocations<MODEL>    ModelAtLocations_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncrement_;
  typedef ObsOperator<MODEL>         ObsOperator_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef Variables<MODEL>           Variables_;

 public:
  static const std::string classname() {return "oops::LinearObsOperator";}

  explicit LinearObsOperator(const ObsOperator_ &);
  explicit LinearObsOperator(const LinearObsOperator &);
  ~LinearObsOperator();

/// Interfacing
  const LinearObsOperator_ & linearobsoperator() const {return *oper_;}

/// Obs Operators
  void setTrajectory(const ModelAtLocations_ &, const ObsAuxControl_ &);
  void obsEquivTL(const ModelAtLocations_ &, ObsVector_ &, const ObsAuxIncrement_ &) const;
  void obsEquivAD(ModelAtLocations_ &, const ObsVector_ &, ObsAuxIncrement_ &) const;

/// Other
  Variables_ variables() const;  // Required inputs variables from LinearModel

 private:
  LinearObsOperator & operator=(const LinearObsOperator &);
  void print(std::ostream &) const;
  boost::shared_ptr<LinearObsOperator_> oper_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperator<MODEL>::LinearObsOperator(const ObsOperator_ & hop): oper_()
{
  Log::trace() << "LinearObsOperator<MODEL>::LinearObsOperator starting" << std::endl;
  util::Timer timer(classname(), "LinearObsOperator");
  oper_.reset(LinearObsOperator_::create(hop.obsoperator()));
  Log::trace() << "LinearObsOperator<MODEL>::LinearObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
LinearObsOperator<MODEL>::LinearObsOperator(const LinearObsOperator & other) : oper_(other.oper_)
{
  Log::trace() << "LinearObsOperator<MODEL>::LinearObsOperator copied" << std::endl;
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
void LinearObsOperator<MODEL>::setTrajectory(const ModelAtLocations_ & gom, const ObsAuxControl_ & aux) {
  Log::trace() << "LinearObsOperator<MODEL>::obsEquiv starting" << std::endl;
  util::Timer timer(classname(), "ObsEquiv");
  oper_->setTrajectory(gom.modelatlocations(), aux.obsauxcontrol());
  Log::trace() << "LinearObsOperator<MODEL>::obsEquiv done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LinearObsOperator<MODEL>::obsEquivTL(const ModelAtLocations_ & gom, ObsVector_ & yy,
                                          const ObsAuxIncrement_ & aux) const {
  Log::trace() << "LinearObsOperator<MODEL>::obsEquivTL starting" << std::endl;
  util::Timer timer(classname(), "ObsEquivTL");
  oper_->obsEquivTL(gom.modelatlocations(), yy.obsvector(), aux.obsauxincrement());
  Log::trace() << "LinearObsOperator<MODEL>::obsEquivTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void LinearObsOperator<MODEL>::obsEquivAD(ModelAtLocations_ & gom, const ObsVector_ & yy,
                                          ObsAuxIncrement_ & aux) const {
  Log::trace() << "LinearObsOperator<MODEL>::obsEquivAD starting" << std::endl;
  util::Timer timer(classname(), "ObsEquivAD");
  oper_->obsEquivAD(gom.modelatlocations(), yy.obsvector(), aux.obsauxincrement());
  Log::trace() << "LinearObsOperator<MODEL>::obsEquivAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL> LinearObsOperator<MODEL>::variables() const {
  Log::trace() << "LinearObsOperator<MODEL>::variables starting" << std::endl;
  util::Timer timer(classname(), "variables");
  Variables_ inputs(oper_->variables());
  Log::trace() << "LinearObsOperator<MODEL>::variables done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LinearObsOperator<MODEL>::print(std::ostream & os) const {
  Log::trace() << "LinearObsOperator<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
//  os << *increment_;
  Log::trace() << "LinearObsOperator<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LINEAROBSOPERATOR_H_
