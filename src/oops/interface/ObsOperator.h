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

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsOperator : public util::Printable,
                    private boost::noncopyable,
                    private util::ObjectCounter<ObsOperator<MODEL> > {
  typedef typename MODEL::ObsOperator  ObsOperator_;
  typedef GeoVaLs<MODEL>               GeoVaLs_;
  typedef Locations<MODEL>             Locations_;
  typedef ObsAuxControl<MODEL>         ObsAuxControl_;
  typedef ObsVector<MODEL>             ObsVector_;
  typedef ObservationSpace<MODEL>      ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObsOperator";}

  ObsOperator(const ObsSpace_ &, const eckit::Configuration &);
  ~ObsOperator();

/// Obs Operator
  void simulateObs(const GeoVaLs_ &, ObsVector_ &, const ObsAuxControl_ &) const;

/// Interfacing
  const ObsOperator_ & obsoperator() const {return *oper_;}

/// Other
  const Variables & variables() const;  // Required input variables from Model
  Locations_ locations(const util::DateTime &, const util::DateTime &) const;
  const std::string & obstype() const {return oper_->obstype();}

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<ObsOperator_> oper_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsOperator<MODEL>::ObsOperator(const ObsSpace_ & os,
                                const eckit::Configuration & config) : oper_() {
  Log::trace() << "ObsOperator<MODEL>::ObsOperator starting" << std::endl;
  util::Timer timer(classname(), "ObsOperator");
  oper_.reset(new ObsOperator_(os.observationspace(), config));
  Log::trace() << "ObsOperator<MODEL>::ObsOperator done" << std::endl;
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
void ObsOperator<MODEL>::simulateObs(const GeoVaLs_ & gvals, ObsVector_ & yy,
                                     const ObsAuxControl_ & aux) const {
  Log::trace() << "ObsOperator<MODEL>::simulateObs starting" << std::endl;
  util::Timer timer(classname(), "simulateObs");
  oper_->simulateObs(gvals.geovals(), yy.obsvector(), aux.obsauxcontrol());
  Log::trace() << "ObsOperator<MODEL>::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
const Variables & ObsOperator<MODEL>::variables() const {
  Log::trace() << "ObsOperator<MODEL>::variables starting" << std::endl;
  util::Timer timer(classname(), "variables");
  return oper_->variables();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Locations<MODEL> ObsOperator<MODEL>::locations(const util::DateTime & t1,
                                               const util::DateTime & t2) const {
  Log::trace() << "ObsOperator<MODEL>::locations starting" << std::endl;
  util::Timer timer(classname(), "locations");
  return Locations_(oper_->locations(t1, t2));
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsOperator<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsOperator<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *oper_;
  Log::trace() << "ObsOperator<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSOPERATOR_H_
