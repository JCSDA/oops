/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_GOMQG_H_
#define QG_MODEL_GOMQG_H_

#include <ostream>
#include <string>

#include "model/QgFortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace oops {
  class Variables;
}

namespace qg {
  class LocationsQG;

/// GomQG class to handle local model values for QG model.

class GomQG : public util::Printable,
              private util::ObjectCounter<GomQG> {
 public:
  static const std::string classname() {return "qg::GomQG";}

  GomQG(const LocationsQG &, const oops::Variables &);
  GomQG(const eckit::Configuration &, const oops::Variables &);

  explicit GomQG(): keyGom_(0) {}
  explicit GomQG(int & fgom): keyGom_(fgom) {}

  ~GomQG();

  void zero();
  void random();
  GomQG & operator*=(const double &);
  double dot_product_with(const GomQG &) const;
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

  int & toFortran() {return keyGom_;}
  const int & toFortran() const {return keyGom_;}

 private:
  void print(std::ostream &) const;
  F90goms keyGom_;
};

}  // namespace qg

#endif  // QG_MODEL_GOMQG_H_
