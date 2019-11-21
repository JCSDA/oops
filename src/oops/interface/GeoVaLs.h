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

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
template <typename MODEL>
class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs<MODEL> > {
  typedef typename MODEL::GeoVaLs          GeoVaLs_;
  typedef ObsSpace<MODEL>                  ObsSpace_;
  typedef Locations<MODEL>                 Locations_;

 public:
  static const std::string classname() {return "oops::GeoVaLs";}

  GeoVaLs(const Locations_ &, const Variables &);
  GeoVaLs(const eckit::Configuration &, const ObsSpace_ &,
          const Variables &);
  GeoVaLs(const GeoVaLs &);

  ~GeoVaLs();

/// Interfacing
  const GeoVaLs_ & geovals() const {return *gvals_;}
  GeoVaLs_ & geovals() {return *gvals_;}

/// Linear algebra and utilities, mostly for writing tests
  void abs();
  void zero();
  void random();
  double norm() const;
  GeoVaLs & operator=(const GeoVaLs &);
  GeoVaLs & operator*=(const double &);
  GeoVaLs & operator+=(const GeoVaLs &);
  GeoVaLs & operator-=(const GeoVaLs &);
  GeoVaLs & operator*=(const GeoVaLs &);
  GeoVaLs & operator/=(const GeoVaLs &);
  double dot_product_with(const GeoVaLs &) const;
  void read(const eckit::Configuration &);
  void analytic_init(const Locations_ &, const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<GeoVaLs_> gvals_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL>::GeoVaLs(const Locations_ & locs, const Variables & vars) : gvals_() {
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(locs.locations(), vars));
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
  GeoVaLs<MODEL>::GeoVaLs(const eckit::Configuration & conf,
                          const ObsSpace_ & ospace, const Variables & vars)
  : gvals_() {
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs read starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(conf, ospace.obsspace(), vars));
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL>::GeoVaLs(const GeoVaLs & other): gvals_() {
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(*other.gvals_));
  Log::trace() << "ObsVector<MODEL>::GeoVaLs done" << std::endl;
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

template <typename MODEL>
GeoVaLs<MODEL> & GeoVaLs<MODEL>::operator=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *gvals_ = *rhs.gvals_;
  Log::trace() << "GeovaLs<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL> & GeoVaLs<MODEL>::operator+=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<MODEL>::+=(GeoVaLs, GeoVaLs) starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *gvals_ += *rhs.gvals_;
  Log::trace() << "GeoVaLs<MODEL>::+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL> & GeoVaLs<MODEL>::operator-=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<MODEL>::-=(GeoVaLs, GeoVaLs) starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *gvals_ -= *rhs.gvals_;
  Log::trace() << "GeoVaLs<MODEL>::-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeoVaLs<MODEL> & GeoVaLs<MODEL>::operator*=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<MODEL>::*=(GeoVaLs, GeoVaLs) starting" << std::endl;
  util::Timer timer(classname(), "operator*=(schur)");
  *gvals_ *= *rhs.gvals_;
  Log::trace() << "GeoVaLs<MODEL>::*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------
/*! GeoVaLs Normalization Operator
 *
 * This is a normalization operator that first computes the normalization
 * factor for each variable based on the rms amplitude of that variable across
 * all locations in the reference GeoVaLs object (rhs).  Then each element of
 * the input GeoVals object (*this) is divided by these normalization factors.
 */

template <typename MODEL>
GeoVaLs<MODEL> & GeoVaLs<MODEL>::operator/=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<MODEL>::/=(GeoVaLs, GeoVaLs) starting" << std::endl;
  util::Timer timer(classname(), "operator/=");
  *gvals_ /= *rhs.gvals_;
  Log::trace() << "GeoVaLs<MODEL>::/= done" << std::endl;
  return *this;
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
void GeoVaLs<MODEL>::abs() {
  Log::trace() << "GeoVaLs<MODEL>::abs starting" << std::endl;
  util::Timer timer(classname(), "abs");
  gvals_->abs();
  Log::trace() << "GeoVaLs<MODEL>::abs done" << std::endl;
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
double GeoVaLs<MODEL>::norm() const {
  Log::trace() << "GeoVaLs<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = gvals_->norm();
  Log::trace() << "GeoVaLs<MODEL>::norm done" << std::endl;
  return zz;
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
/*! \brief GeoVaLs Analytic Initialization
 *
 * \details **analytic_init()** was introduced in May, 2018 (initially as a
 * constructor) for use with the interpolation test.  The interpolation test
 * requires an initialization of a GeoVaLs object based on the same analytic
 * formulae used for the State initialization (see test::TestStateInterpolation()
 * for further information).  This in turn requires information about the
 * vertical profile in addition to the latitude and longitude positional
 * information in the Locations object.  Currently, this information
 * about the vertical profile is obtained from an existing GeoVaLs object
 * (passed as *other*) that represents the output of the State::interpolate()
 * method.  The State.StateGenerate section of the configuration file is
 * also passed to this constructor to provide further information required
 * for the analytic initialization.
 *
 * \date May, 2018: created as a constructor (M. Miesch, JCSDA)
 * \date June, 2018: moved to a method (M. Miesch, JCSDA)
 *
 * \sa test::TestStateInterpolation()
 */

template <typename MODEL>
  void GeoVaLs<MODEL>::analytic_init(const Locations_ & locs,
                          const eckit::Configuration & conf) {
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs analytic init starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_->analytic_init(locs.locations(), conf);
  Log::trace() << "GeoVaLs<MODEL>::GeoVaLs analytic init done" << std::endl;
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
