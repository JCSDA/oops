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

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/GeometryQGIterator.h"
#include "oops/qg/QgFortran.h"

namespace eckit {
  class Configuration;
}

namespace qg {
class GeometryQGIterator;

// -----------------------------------------------------------------------------
/// GeometryQG handles geometry for QG model.

class GeometryQG : public util::Printable,
                   private util::ObjectCounter<GeometryQG> {
 public:
  static const std::string classname() {return "qg::GeometryQG";}

  GeometryQG(const eckit::Configuration &, const eckit::mpi::Comm &);
  GeometryQG(const GeometryQG &);
  ~GeometryQG();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}

  GeometryQGIterator begin() const;
  GeometryQGIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}

 private:
  GeometryQG & operator=(const GeometryQG &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::StructuredColumns> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_GEOMETRYQG_H_
