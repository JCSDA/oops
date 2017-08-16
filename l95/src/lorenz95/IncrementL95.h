/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_INCREMENTL95_H_
#define LORENZ95_INCREMENTL95_H_

#include <ostream>
#include <string>
#include <vector>

#include "lorenz95/FieldL95.h"
#include "lorenz95/Resolution.h"

#include "oops/base/GeneralizedDepartures.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
}

namespace lorenz95 {
  class GomL95;
  class LocsL95;
  class ModelBiasCorrection;
  class StateL95;
  class NoVariables;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------
class IncrementL95 : public util::Printable,
                     public oops::GeneralizedDepartures,
                     private util::ObjectCounter<IncrementL95> {
 public:
  static const std::string classname() {return "lorenz95::IncrementL95";}

/// Constructor, destructor
  IncrementL95(const Resolution &, const NoVariables &, const util::DateTime &);
  IncrementL95(const Resolution &, const IncrementL95 &);
  IncrementL95(const IncrementL95 &, const bool);
  virtual ~IncrementL95();

/// Basic operators
  void diff(const StateL95 &, const StateL95 &);
  void zero();
  void zero(const util::DateTime &);
  void dirac(const eckit::Configuration &);
  IncrementL95 & operator =(const IncrementL95 &);
  IncrementL95 & operator+=(const IncrementL95 &);
  IncrementL95 & operator-=(const IncrementL95 &);
  IncrementL95 & operator*=(const double &);
  void axpy(const double &, const IncrementL95 &, const bool check = true);
  double dot_product_with(const IncrementL95 &) const;
  void schur_product_with(const IncrementL95 &);
  void random();

/// Interpolate to observation location
  void interpolateTL(const LocsL95 &, GomL95 &) const;
  void interpolateAD(const LocsL95 &, const GomL95 &);

/// Convert to/from generic unstructured grid
  void convert_to(oops::UnstructuredGrid &) const;
  void convert_from(const oops::UnstructuredGrid &);

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm () const {return fld_.rms();}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

/// Access to data
  const FieldL95 & getField() const {return fld_;}
  FieldL95 & getField() {return fld_;}
  boost::shared_ptr<const Resolution> geometry() const {
    boost::shared_ptr<const Resolution> geom(new Resolution(fld_.resol()));
    return geom;
  }
  std::vector<double> & asVector() {return fld_.asVector();}
  const std::vector<double> & asVector() const {return fld_.asVector();}

  void accumul(const double &, const StateL95 &);

 private:
  void print(std::ostream &) const;
  FieldL95 fld_;
  util::DateTime time_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_INCREMENTL95_H_
