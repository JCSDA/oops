/*
 * (C) Copyright 2018  UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_INTERFACE_OBSDIAGNOSTICS_H_
#define OOPS_INTERFACE_OBSDIAGNOSTICS_H_

#include <ostream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsDiagnostics : public util::Printable,
                       private boost::noncopyable,
                       private util::ObjectCounter<ObsDiagnostics<MODEL> > {
  typedef typename MODEL::ObsDiagnostics   ObsDiags_;
  typedef ObservationSpace<MODEL>          ObsSpace_;
  typedef Locations<MODEL>                 Locations_;

 public:
  static const std::string classname() {return "oops::ObsDiagnostics";}

  ObsDiagnostics(const ObsSpace_ &, const Locations_ &, const Variables &);
  ObsDiagnostics(const eckit::Configuration &, const ObsSpace_ &, const Variables &);

  ~ObsDiagnostics();

/// Interfacing
  ObsDiags_ & obsdiagnostics() {return *diags_;}
  const ObsDiags_ & obsdiagnostics() const {return *diags_;}

// I/O
  void save(const std::string &) const;

 private:
  void print(std::ostream &) const;
  boost::shared_ptr<ObsDiags_> diags_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsDiagnostics<MODEL>::ObsDiagnostics(const ObsSpace_ & os, const Locations_ & locs,
                                      const Variables & vars) : diags_()
{
  Log::trace() << "ObsDiagnostics<MODEL>::ObsDiagnostics starting" << std::endl;
  util::Timer timer(classname(), "ObsDiagnostics");
  diags_.reset(new ObsDiags_(os.observationspace(), locs.locations(), vars));
  Log::trace() << "ObsDiagnostics<MODEL>::ObsDiagnostics done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsDiagnostics<MODEL>::ObsDiagnostics(const eckit::Configuration & conf, const ObsSpace_ & os,
                                      const Variables & vars) : diags_()
{
  Log::trace() << "ObsDiagnostics<MODEL>::ObsDiagnostics starting" << std::endl;
  util::Timer timer(classname(), "ObsDiagnostics");
  diags_.reset(new ObsDiags_(conf, os.observationspace(), vars));
  Log::trace() << "ObsDiagnostics<MODEL>::ObsDiagnostics done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
ObsDiagnostics<MODEL>::~ObsDiagnostics() {
  Log::trace() << "ObsDiagnostics<MODEL>::~ObsDiagnostics starting" << std::endl;
  util::Timer timer(classname(), "~ObsDiagnostics");
  diags_.reset();
  Log::trace() << "ObsDiagnostics<MODEL>::~ObsDiagnostics done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsDiagnostics<MODEL>::save(const std::string & name) const {
  Log::trace() << "ObsDiagnostics<MODEL, DATATYPE>::save starting " << name << std::endl;
  util::Timer timer(classname(), "save");
  diags_->save(name);
  Log::trace() << "ObsDiagnostics<MODEL, DATATYPE>::save done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsDiagnostics<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsDiagnostics<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *diags_;
  Log::trace() << "ObsDiagnostics<MODEL>::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSDIAGNOSTICS_H_
