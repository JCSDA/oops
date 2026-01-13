/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_GOML95_H_
#define LORENZ95_GOML95_H_

#include <Eigen/Core>

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variable;
  class Variables;
  template <typename OBS> class Locations;
}

namespace lorenz95 {
  struct L95ObsTraits;
  class ObsTable;


/// GomL95 class to handle State values at obs locations for L95 model.
class GomL95 : public util::Printable,
               private util::ObjectCounter<GomL95> {
  /// References to read-only or writable vector- or matrix-valued expressions.
  ///
  /// For example, an Eigen::Vector, Eigen::Matrix or an Eigen::Map (the latter can be used as a
  /// view onto a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using ConstVectorRef = Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>;
  template <typename T>
  using ConstMatrixRef = Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;
  template <typename T>
  using MatrixRef = Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

 public:
  static const std::string classname() {return "lorenz95::GomL95";}

  GomL95(const oops::Locations<L95ObsTraits> & locs,
         const oops::Variables & vars, const std::vector<size_t> & sizes);
  GomL95(const eckit::Configuration &, const ObsTable &, const oops::Variables &);

  void zero();
  void random();
  double rms() const;
  double normalizedrms(const GomL95 &) const;
  GomL95 & operator*=(const double &);
  GomL95 & operator+=(const GomL95 &);
  GomL95 & operator-=(const GomL95 &);
  GomL95 & operator*=(const GomL95 &);
  double dot_product_with(const GomL95 &) const;
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  void print(std::ostream &) const;

  size_t size() const {return size_;}
  const double & operator[](const int ii) const {return locval_[ii];}
  double & operator[](const int ii) {return locval_[ii];}

  void fill(const oops::Variable &name, const ConstVectorRef<size_t> &indx,
            const ConstMatrixRef<double> &vals, const bool levelsTopDown);
  void fillAD(const oops::Variable &name, const ConstVectorRef<size_t> &indx,
              MatrixRef<double> vals, const bool levelsTopDow) const;

 private:
  size_t size_;
  std::vector<double> locval_;
};

}  // namespace lorenz95

#endif  // LORENZ95_GOML95_H_
