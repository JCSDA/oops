/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSOPERATOR_H_
#define OOPS_INTERFACE_OBSOPERATOR_H_

#include <string>

#include <boost/shared_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObservationSpace.h"
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
class ObsOperator : public util::Printable,
                    private util::ObjectCounter<ObsOperator<MODEL> > {
  typedef typename MODEL::ObsOperator           ObsOperator_;
  typedef ModelAtLocations<MODEL>    ModelAtLocations_;
  typedef ObsAuxControl<MODEL>       ObsAuxControl_;
  typedef ObsVector<MODEL>           ObsVector_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef Variables<MODEL>           Variables_;

 public:
  static const std::string classname() {return "oops::ObsOperator";}

  ObsOperator(const ObsSpace_ &, const eckit::Configuration &);
  ObsOperator(const ObsOperator &);
  ~ObsOperator();

/// Interfacing
  const ObsOperator_ & obsoperator() const {return *oper_;}

/// Obs Operator
  void obsEquiv(const ModelAtLocations_ &, ObsVector_ &, const ObsAuxControl_ &) const;

/// Other
  Variables_ variables() const;  // Required inputs variables from Model
  void generateObsError(const eckit::Configuration &);

 private:
  ObsOperator & operator=(const ObsOperator &);
  void print(std::ostream &) const;
  boost::shared_ptr<ObsOperator_> oper_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperator<MODEL>::ObsOperator(const ObsSpace_ & os, const eckit::Configuration & conf)
  : oper_()
{
  Log::trace() << "ObsOperator<MODEL>::ObsOperator starting" << std::endl;
  util::Timer timer(classname(), "ObsOperator");
  oper_.reset(ObsOperator_::create(os.observationspace(), conf));
  Log::trace() << "ObsOperator<MODEL>::ObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperator<MODEL>::ObsOperator(const ObsOperator & other) : oper_(other.oper_)
{
  Log::trace() << "ObsOperator<MODEL>::ObsOperator copied" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperator<MODEL>::~ObsOperator() {
  Log::trace() << "ObsOperator<MODEL>::~ObsOperator starting" << std::endl;
  util::Timer timer(classname(), "~ObsOperator");
  oper_.reset();
  Log::trace() << "ObsOperator<MODEL>::~ObsOperator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsOperator<MODEL>::obsEquiv(const ModelAtLocations_ & gom, ObsVector_ & yy,
                                  const ObsAuxControl_ & aux) const {
  Log::trace() << "ObsOperator<MODEL>::obsEquiv starting" << std::endl;
  util::Timer timer(classname(), "ObsEquiv");
  oper_->obsEquiv(gom.modelatlocations(), yy.obsvector(), aux.obsauxcontrol());
  Log::trace() << "ObsOperator<MODEL>::obsEquiv done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Variables<MODEL> ObsOperator<MODEL>::variables() const {
  Log::trace() << "ObsOperator<MODEL>::variables starting" << std::endl;
  util::Timer timer(classname(), "variables");
  Variables<MODEL> var(oper_->variables());
  Log::trace() << "ObsOperator<MODEL>::variables done" << std::endl;
  return var;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsOperator<MODEL>::generateObsError(const eckit::Configuration & conf) {
  Log::trace() << "ObsOperator<MODEL>::generateObsError starting" << std::endl;
  util::Timer timer(classname(), "generateObsError");
  oper_->generateObsError(conf);
  Log::trace() << "ObsOperator<MODEL>::generateObsError done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsOperator<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsOperator<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
//  os << *increment_;
  Log::trace() << "ObsOperator<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSOPERATOR_H_
