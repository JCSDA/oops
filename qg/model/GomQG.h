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

#include <Eigen/Core>

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "oops/qg/QgFortran.h"

namespace oops {
  template <typename OBS> class Locations;
  class Variables;
}

namespace qg {
  struct QgObsTraits;

/// \brief Parameters controlling a QG GeoVaLs read/write
class GomQGParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GomQGParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> filename{"filename", "filename for input and output",
                                                this};
};

/// GomQG class to handle local model values for QG model.

class GomQG : public util::Printable,
              private util::ObjectCounter<GomQG> {
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
  typedef GomQGParameters Parameters_;
  typedef oops::Locations<QgObsTraits> Locations_;

  static const std::string classname() {return "qg::GomQG";}

  GomQG(const Locations_ & locs,
        const oops::Variables & vars, const std::vector<size_t> & sizes);
  GomQG(const Parameters_ & params, const ObsSpaceQG & ospace, const oops::Variables & vars);
  explicit GomQG(const GomQG &);

  GomQG(): keyGom_(0) {}
  explicit GomQG(int & fgom): keyGom_(fgom) {}

  ~GomQG();

  void zero();
  void random();
  double rms() const;
  double normalizedrms(const GomQG &) const;
  GomQG & operator=(const GomQG &);
  GomQG & operator*=(const double &);
  GomQG & operator+=(const GomQG &);
  GomQG & operator-=(const GomQG &);
  GomQG & operator*=(const GomQG &);
  double dot_product_with(const GomQG &) const;
  void read(const Parameters_ &);
  void write(const Parameters_ &) const;

  const int & toFortran() const {return keyGom_;}

  void fill(const std::string &name, const ConstVectorRef<size_t> &indx,
            const ConstMatrixRef<double> &vals, const bool levelsTopDown);
  void fillAD(const std::string &name, const ConstVectorRef<size_t> &indx,
              MatrixRef<double> vals, const bool levelsTopDown) const;

 private:
  void print(std::ostream &) const;
  F90gom keyGom_;
  oops::Variables vars_;
};

}  // namespace qg

#endif  // QG_MODEL_GOMQG_H_
