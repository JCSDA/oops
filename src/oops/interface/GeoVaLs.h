/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_GEOVALS_H_
#define OOPS_INTERFACE_GEOVALS_H_

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class GeoVaLs : public util::Printable,
                private boost::noncopyable,
                private util::ObjectCounter<GeoVaLs<MODEL> > {
  typedef typename MODEL::GeoVaLs          GeoVaLs_;
  typedef Geometry<MODEL>                  Geometry_;
  typedef ObservationSpace<MODEL>          ObsSpace_;
  typedef Variables<MODEL>                 Variables_;

 public:
  static const std::string classname() {return "oops::GeoVaLs";}

  GeoVaLs(const ObsSpace_ &, const Variables_ &,
          const util::DateTime &, const util::DateTime &,
          const Geometry_ &);
  explicit GeoVaLs(const eckit::Configuration &);
  ~GeoVaLs();

/// Interfacing
  const GeoVaLs_ & geovals() const {return *gvals_;}
  GeoVaLs_ & geovals() {return *gvals_;}

/// Linear algebra and utilities, mostly for writing tests
  void zero();
  void random();
  GeoVaLs & operator*=(const double &);
  double dot_product_with(const GeoVaLs &) const;
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<GeoVaLs_> gvals_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL>::GeoVaLs(const ObsSpace_ & os, const Variables_ & var,
                        const util::DateTime & t1, const util::DateTime & t2,
                        const Geometry_ & resol) : gvals_() {
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(os.observationspace(), var.variables(), t1, t2, resol.geometry()));
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL>::GeoVaLs(const eckit::Configuration & conf) : gvals_() {
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs read starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(conf));
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL>::~GeoVaLs() {
  Log::trace() << "GeoVaLs<MODEL>::~GeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "~GeoVaLs");
  gvals_.reset();
  Log::trace() << "GeoVaLs<MODEL>::~GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
double GeoVaLs<MODEL>::dot_product_with(const GeoVaLs & other) const {
  Log::trace() << "GeoVaLs<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = gvals_->dot_product_with(*other.gvals_);
  Log::trace() << "GeoVaLs<MODEL>::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
GeoVaLs<MODEL> & GeoVaLs<MODEL>::operator*=(const double & zz) {
  Log::trace() << "GeoVaLs<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *gvals_ *= zz;
  Log::trace() << "GeoVaLs<MODEL>::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void GeoVaLs<MODEL>::zero() {
  Log::trace() << "GeoVaLs<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  gvals_->zero();
  Log::trace() << "GeoVaLs<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void GeoVaLs<MODEL>::random() {
  Log::trace() << "GeoVaLs<MODEL>::random starting" << std::endl;
  util::Timer timer(classname(), "random");
  gvals_->random();
  Log::trace() << "GeoVaLs<MODEL>::random done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GeoVaLs<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "GeoVaLs<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  gvals_->read(conf);
  Log::trace() << "GeoVaLs<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GeoVaLs<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "GeoVaLs<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  gvals_->write(conf);
  Log::trace() << "GeoVaLs<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GeoVaLs<MODEL>::print(std::ostream & os) const {
  Log::trace() << "GeoVaLs<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *gvals_;
  Log::trace() << "GeoVaLs<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_GEOVALS_H_
