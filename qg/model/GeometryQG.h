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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "oops/qg/GeometryQGIterator.h"
#include "oops/qg/QgFortran.h"

namespace oops {
  class Variables;
}

namespace qg {

class GeometryQgParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryQgParameters, Parameters)

 public:
  /// Domain size
  oops::RequiredParameter<int> nx{"nx", this};
  oops::RequiredParameter<int> ny{"ny", this};
  /// Depths
  oops::RequiredParameter<std::vector<float>> depths{"depths", this};
  /// Heating option (AS: should it be in geometry or model?)
  oops::Parameter<bool> heating{"heating", true, this};
};

class GeometryQGIterator;

// -----------------------------------------------------------------------------
/// GeometryQG handles geometry for QG model.

class GeometryQG : public util::Printable,
                   private util::ObjectCounter<GeometryQG> {
 public:
  typedef GeometryQgParameters Parameters_;

  static const std::string classname() {return "qg::GeometryQG";}

  GeometryQG(const GeometryQgParameters &, const eckit::mpi::Comm &);
  GeometryQG(const GeometryQG &);
  ~GeometryQG();

  const F90geom & toFortran() const {return keyGeom_;}

  GeometryQGIterator begin() const;
  GeometryQGIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}
  size_t levels() const {return levs_;}

  std::vector<size_t> variableSizes(const oops::Variables & vars) const;

  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;

 private:
  GeometryQG & operator=(const GeometryQG &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
  size_t levs_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_GEOMETRYQG_H_
