/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_UNSTRUCTUREDGRID_H_
#define OOPS_GENERIC_UNSTRUCTUREDGRID_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

class UnstructuredGrid : public util::Printable,
                         private boost::noncopyable,
                         private util::ObjectCounter<UnstructuredGrid> {
 public:
  static const std::string classname() {return "oops::UnstructuredGrid";}

  explicit UnstructuredGrid(const int &, const int &);
  explicit UnstructuredGrid(UnstructuredGrid &);
  ~UnstructuredGrid();

// Will be useful for tests
  void zero();
  void random();
  double dot_product_with(const UnstructuredGrid &) const;

// Will be useful for tests
  int & toFortran() {return keyUGrid_;}
  const int & toFortran() const {return keyUGrid_;}

/// ATLAS-related methods
  void defineGeometry();
  void defineGrids(std::vector<eckit::LocalConfiguration> &) const;
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

 private:
  void print(std::ostream &) const;

  int keyUGrid_;
  std::unique_ptr<atlas::FunctionSpace> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_UNSTRUCTUREDGRID_H_
