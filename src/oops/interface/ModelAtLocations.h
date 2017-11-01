/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_MODELATLOCATIONS_H_
#define OOPS_INTERFACE_MODELATLOCATIONS_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/interface/Geometry.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/Variables.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace util {
  class DateTime;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ModelAtLocations : public util::Printable,
                         private boost::noncopyable,
                         private util::ObjectCounter<ModelAtLocations<MODEL> > {
  typedef typename MODEL::ModelAtLocations GOM_;
  typedef Geometry<MODEL>                  Geometry_;
  typedef ObservationSpace<MODEL>          ObsSpace_;
  typedef Variables<MODEL>                 Variables_;

 public:
  static const std::string classname() {return "oops::ModelAtLocations";}

// Geometry passed to constructor for IFS... :-(
  ModelAtLocations(const ObsSpace_ &, const Variables_ &,
                   const util::DateTime &, const util::DateTime &,
                   const Geometry_ &);
  ~ModelAtLocations();

  double dot_product_with(const ModelAtLocations &) const;
  void zero();

/// Interfacing
  const GOM_ & modelatlocations() const {return *gom_;}
  GOM_ & modelatlocations() {return *gom_;}

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<GOM_> gom_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelAtLocations<MODEL>::ModelAtLocations(const ObsSpace_ & os, const Variables_ & var,
                                          const util::DateTime & t1, const util::DateTime & t2,
                                          const Geometry_ & resol)
  : gom_()
{
  Log::trace() << "ModelAtLocations<MODEL>::ModelAtLocations starting" << std::endl;
  util::Timer timer(classname(), "ModelAtLocations");
  gom_.reset(new GOM_(os.observationspace(), var.variables(), t1, t2, resol.geometry()));
  Log::trace() << "ModelAtLocations<MODEL>::ModelAtLocations done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelAtLocations<MODEL>::~ModelAtLocations() {
  Log::trace() << "ModelAtLocations<MODEL>::~ModelAtLocations starting" << std::endl;
  util::Timer timer(classname(), "~ModelAtLocations");
  gom_.reset();
  Log::trace() << "ModelAtLocations<MODEL>::~ModelAtLocations done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
double ModelAtLocations<MODEL>::dot_product_with(const ModelAtLocations & other) const {
  Log::trace() << "ModelAtLocations<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = gom_->dot_product_with(other.gom_);
  Log::trace() << "ModelAtLocations<MODEL>::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ModelAtLocations<MODEL>::zero() {
  Log::trace() << "ModelAtLocations<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  gom_->zero();
  Log::trace() << "ModelAtLocations<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAtLocations<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ModelAtLocations<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
//  os << *increment_;
  Log::trace() << "ModelAtLocations<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELATLOCATIONS_H_
