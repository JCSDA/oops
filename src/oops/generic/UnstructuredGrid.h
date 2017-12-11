/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_UNSTRUCTURED_GRID_H_
#define OOPS_GENERIC_UNSTRUCTURED_GRID_H_

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

class UnstructuredGrid : public util::Printable,
                         private boost::noncopyable,
                         private util::ObjectCounter<UnstructuredGrid> {
 public:
  static const std::string classname() {return "oops::UnstructuredGrid";}

  UnstructuredGrid();
  ~UnstructuredGrid();

// Get local geometry
  std::vector<double> getLats();
  std::vector<double> getLons();
  std::vector<double> getAreas();
  std::vector<double> getVunit();
  std::vector<int> getMask(const int &);
  std::vector<int> getGlbInd();
  int getNvar();
  std::vector<double> getData();

// Will be useful for tests
  void zero();
  void random();
  double dot_product_with(const UnstructuredGrid &) const;

// Will be useful for tests
  int & toFortran() {return keyUGrid_;}
  const int & toFortran() const {return keyUGrid_;}

 private:
  void print(std::ostream &) const;

  int keyUGrid_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_UNSTRUCTURED_GRID_H_
