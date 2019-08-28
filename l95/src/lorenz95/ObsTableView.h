/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSTABLEVIEW_H_
#define LORENZ95_OBSTABLEVIEW_H_

#include <ostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/geometry/Point2.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsTable.h"

namespace lorenz95 {

class ObsVec1D;

/// A Simple Observation Data Handler
/*!
 *  ObsTableView defines a simple observation handler
 *  that mimicks the interfaces required from ODB.
 */

// -----------------------------------------------------------------------------
class ObsTableView : public util::Printable,
                     private util::ObjectCounter<ObsTableView> {
 public:
  static const std::string classname() {return "lorenz95::ObsTableView";}

  ObsTableView(const eckit::Configuration &, const util::DateTime &, const util::DateTime &);
  ObsTableView(const ObsTableView &, const eckit::geometry::Point2 &,
               const double &, const int &);
  ~ObsTableView();

  bool has(const std::string &) const;
  void putdb(const std::string &, const std::vector<int> &) const;
  void putdb(const std::string &, const std::vector<float> &) const;
  void putdb(const std::string &, const std::vector<double> &) const;
  void getdb(const std::string &, std::vector<int> &) const;
  void getdb(const std::string &, std::vector<float> &) const;
  void getdb(const std::string &, std::vector<double> &) const;

  void random(std::vector<double> &) const;
  unsigned int nobs() const;
  std::vector<double> locations() const;
  void generateDistribution(const eckit::Configuration &);
  LocsL95 * locations(const util::DateTime & t1, const util::DateTime & t2) const;
  void printJo(const ObsVec1D &, const ObsVec1D &);

  const std::string & obsname() const {return obstable_->obsname();}

  const util::DateTime & windowStart() const {return obstable_->windowStart();}
  const util::DateTime & windowEnd() const {return obstable_->windowEnd();}
  const oops::Variables & obsvariables() const { return obstable_->obsvariables(); }
 private:
  void print(std::ostream &) const;
  boost::shared_ptr<ObsTable> obstable_;
  std::vector<int> localobs_;
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSTABLEVIEW_H_
