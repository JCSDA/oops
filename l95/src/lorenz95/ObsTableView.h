/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSTABLEVIEW_H_
#define LORENZ95_OBSTABLEVIEW_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/geometry/Point2.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsTable.h"

namespace lorenz95 {

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

  ObsTableView(const eckit::Configuration &, const eckit::mpi::Comm &,
               const util::DateTime &, const util::DateTime &, const eckit::mpi::Comm &);
  ObsTableView(const ObsTableView &, const eckit::geometry::Point2 &,
               const eckit::Configuration &);
  ~ObsTableView();

  void putdb(const std::string &, const std::vector<int> &) const;
  void putdb(const std::string &, const std::vector<float> &) const;
  void putdb(const std::string &, const std::vector<double> &) const;
  void getdb(const std::string &, std::vector<int> &) const;
  void getdb(const std::string &, std::vector<float> &) const;
  void getdb(const std::string &, std::vector<double> &) const;

  void random(std::vector<double> &) const;
  unsigned int nobs() const;
  void generateDistribution(const eckit::Configuration &);
  std::unique_ptr<LocsL95> locations() const;

  size_t index(const size_t ii) const {return localobs_[ii];}
  const std::string & obsname() const {return obstable_->obsname();}

  const util::DateTime & windowStart() const {return obstable_->windowStart();}
  const util::DateTime & windowEnd() const {return obstable_->windowEnd();}
  const oops::Variables & obsvariables() const { return obstable_->obsvariables(); }
  const std::vector<double> & obsdist() const {return obsdist_;}
 private:
  void print(std::ostream &) const;
  std::shared_ptr<ObsTable> obstable_;
  std::vector<size_t> localobs_;
  std::vector<double> obsdist_;
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSTABLEVIEW_H_
