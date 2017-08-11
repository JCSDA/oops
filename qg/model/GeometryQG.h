/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_GEOMETRYQG_H_
#define QG_MODEL_GEOMETRYQG_H_

#include <ostream>
#include <string>
#include <vector>

#include "model/QgFortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace qg {

// -----------------------------------------------------------------------------
/// GeometryQG handles geometry for QG model.

class GeometryQG : public util::Printable,
                   private util::ObjectCounter<GeometryQG> {
 public:
  static const std::string classname() {return "qg::GeometryQG";}

  explicit GeometryQG(const eckit::Configuration &);
  GeometryQG(const GeometryQG &);
  ~GeometryQG();

  std::vector<int> getDims() const;
  std::vector<double> getLats() const;
  std::vector<double> getLons() const;
  std::vector<double> getLevs() const;
  std::vector<double> getArea() const;
  std::vector<int> getMask(const int &) const;

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}

 private:
  GeometryQG & operator=(const GeometryQG &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_GEOMETRYQG_H_
