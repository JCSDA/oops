/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_GEOVALS_H_
#define OOPS_INTERFACE_GEOVALS_H_

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------
template <typename OBS>
class GeoVaLs : public util::Printable,
                private util::ObjectCounter<GeoVaLs<OBS> > {
  typedef typename OBS::GeoVaLs          GeoVaLs_;
  typedef ObsSpace<OBS>                  ObsSpace_;
  typedef Locations<OBS>                 Locations_;

  /// \brief A reference to a read-only vector-valued expression.
  ///
  /// For example, an Eigen::Vector or an Eigen::Map (the latter can be used as a view onto
  /// a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using ConstVectorRef = Eigen::Ref<const Eigen::Vector<T, Eigen::Dynamic>>;

  /// \brief A reference to a read-only matrix-valued expression.
  ///
  /// For example, an Eigen::Matrix or an Eigen::Map (the latter can be used as a view onto
  /// a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using ConstMatrixRef = Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

  /// \brief A reference to a writable matrix-valued expression.
  ///
  /// For example, an Eigen::Matrix or an Eigen::Map (the latter can be used as a view onto
  /// a chunk of memory stored in another container, such as a std::vector).
  template <typename T>
  using MatrixRef = Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>;

 public:
  typedef typename GeoVaLs_::Parameters_ Parameters_;

  static const std::string classname() {return "oops::GeoVaLs";}

  /// Allocate GeoVaLs for \p locs locations, to be filled with \p vars variables.
  /// Sizes of GeoVaLs for i-th variable at a single location are defined by
  /// i-th value of \p sizes.
  GeoVaLs(const Locations_ & locs, const Variables &, const std::vector<size_t> & sizes);
  GeoVaLs(const Parameters_ &, const ObsSpace_ &, const Variables &);
  GeoVaLs(const GeoVaLs &);

  ~GeoVaLs();

/// Interfacing
  const GeoVaLs_ & geovals() const {return *gvals_;}
  GeoVaLs_ & geovals() {return *gvals_;}

/// Linear algebra and utilities, mostly for writing tests
  void zero();
  void random();
  double rms() const;
  double normalizedrms(const GeoVaLs &) const;
  GeoVaLs & operator=(const GeoVaLs &);
  GeoVaLs & operator*=(const double &);
  GeoVaLs & operator+=(const GeoVaLs &);
  GeoVaLs & operator-=(const GeoVaLs &);
  GeoVaLs & operator*=(const GeoVaLs &);
  double dot_product_with(const GeoVaLs &) const;
  void read(const Parameters_ &);
  void write(const Parameters_ &) const;

  /// \brief Set the values of a given variable at specified locations.
  ///
  /// \param name
  ///   Variable name.
  /// \param indx
  ///   Location indices.
  /// \param vals
  ///   Matrix whose `i`th row contains the values of the variable `name` at the location with
  ///   index `indx[i]`, ordered from top to bottom if `levelsTopDown` is `true` and from bottom to
  ///   top otherwise.
  /// \param levelTopDown
  ///   True if each row of `vals` contains variable values at model levels ordered from top to
  ///   bottom, false if they are ordered from bottom to top.
  void fill(const std::string &name, const ConstVectorRef<size_t> &indx,
            const ConstMatrixRef<double> &vals, const bool levelsTopDown);
  /// \brief Adjoint of fill().
  void fillAD(const std::string &name, const ConstVectorRef<size_t> &indx,
              MatrixRef<double> vals, const bool levelsTopDown) const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<GeoVaLs_> gvals_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
GeoVaLs<OBS>::GeoVaLs(const Locations_ & locs, const Variables & vars,
                      const std::vector<size_t> & sizes) : gvals_() {
  Log::trace() << "GeoVaLs<OBS>::GeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(locs.locations(), vars, sizes));
  Log::trace() << "GeoVaLs<OBS>::GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
  GeoVaLs<OBS>::GeoVaLs(const Parameters_ & params,
                        const ObsSpace_ & ospace, const Variables & vars)
  : gvals_() {
  Log::trace() << "GeoVaLs<OBS>::GeoVaLs read starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(params, ospace.obsspace(), vars));
  Log::trace() << "GeoVaLs<OBS>::GeoVaLs read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
GeoVaLs<OBS>::GeoVaLs(const GeoVaLs & other): gvals_() {
  Log::trace() << "GeoVaLs<OBS>::GeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "GeoVaLs");
  gvals_.reset(new GeoVaLs_(*other.gvals_));
  Log::trace() << "ObsVector<OBS>::GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
GeoVaLs<OBS>::~GeoVaLs() {
  Log::trace() << "GeoVaLs<OBS>::~GeoVaLs starting" << std::endl;
  util::Timer timer(classname(), "~GeoVaLs");
  gvals_.reset();
  Log::trace() << "GeoVaLs<OBS>::~GeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
double GeoVaLs<OBS>::dot_product_with(const GeoVaLs & other) const {
  Log::trace() << "GeoVaLs<OBS>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = gvals_->dot_product_with(*other.gvals_);
  Log::trace() << "GeoVaLs<OBS>::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template <typename OBS>
GeoVaLs<OBS> & GeoVaLs<OBS>::operator=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *gvals_ = *rhs.gvals_;
  Log::trace() << "GeovaLs<OBS>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
GeoVaLs<OBS> & GeoVaLs<OBS>::operator+=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<OBS>::+=(GeoVaLs, GeoVaLs) starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *gvals_ += *rhs.gvals_;
  Log::trace() << "GeoVaLs<OBS>::+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
GeoVaLs<OBS> & GeoVaLs<OBS>::operator-=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<OBS>::-=(GeoVaLs, GeoVaLs) starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *gvals_ -= *rhs.gvals_;
  Log::trace() << "GeoVaLs<OBS>::-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
GeoVaLs<OBS> & GeoVaLs<OBS>::operator*=(const GeoVaLs & rhs) {
  Log::trace() << "GeoVaLs<OBS>::*=(GeoVaLs, GeoVaLs) starting" << std::endl;
  util::Timer timer(classname(), "operator*=(schur)");
  *gvals_ *= *rhs.gvals_;
  Log::trace() << "GeoVaLs<OBS>::*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename OBS>
GeoVaLs<OBS> & GeoVaLs<OBS>::operator*=(const double & zz) {
  Log::trace() << "GeoVaLs<OBS>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *gvals_ *= zz;
  Log::trace() << "GeoVaLs<OBS>::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
double GeoVaLs<OBS>::rms() const {
  Log::trace() << "GeoVaLs<OBS>::rms starting" << std::endl;
  util::Timer timer(classname(), "rms");
  double zz = gvals_->rms();
  Log::trace() << "GeoVaLs<OBS>::rms done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template <typename OBS>
double GeoVaLs<OBS>::normalizedrms(const GeoVaLs & rhs) const {
  Log::trace() << "GeoVaLs<OBS>::normalizedrms starting" << std::endl;
  util::Timer timer(classname(), "normalizedrms");
  double zz = gvals_->normalizedrms(*rhs.gvals_);
  Log::trace() << "GeoVaLs<OBS>::normalizedrms done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void GeoVaLs<OBS>::zero() {
  Log::trace() << "GeoVaLs<OBS>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  gvals_->zero();
  Log::trace() << "GeoVaLs<OBS>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void GeoVaLs<OBS>::random() {
  Log::trace() << "GeoVaLs<OBS>::random starting" << std::endl;
  util::Timer timer(classname(), "random");
  gvals_->random();
  Log::trace() << "GeoVaLs<OBS>::random done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void GeoVaLs<OBS>::fill(const std::string &name, const ConstVectorRef<size_t> &indx,
                        const ConstMatrixRef<double> &vals, const bool levelsTopDown) {
  Log::trace() << "GeoVaLs<OBS>::fill starting" << std::endl;
  util::Timer timer(classname(), "fill");
  ASSERT(indx.size() == vals.rows());
  gvals_->fill(name, indx, vals, levelsTopDown);
  Log::trace() << "GeoVaLs<OBS>::fill done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void GeoVaLs<OBS>::fillAD(const std::string & name,
                          const Eigen::Ref<const Eigen::VectorX<size_t>> &indx,
                          Eigen::Ref<Eigen::MatrixXd> vals,
                          const bool levelsTopDown) const {
  Log::trace() << "GeoVaLs<OBS>::fillAD starting" << std::endl;
  util::Timer timer(classname(), "fillAD");
  ASSERT(indx.size() == vals.rows());
  gvals_->fillAD(name, indx, vals, levelsTopDown);
  Log::trace() << "GeoVaLs<OBS>::fillAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void GeoVaLs<OBS>::read(const Parameters_ & params) {
  Log::trace() << "GeoVaLs<OBS>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  gvals_->read(params);
  Log::trace() << "GeoVaLs<OBS>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void GeoVaLs<OBS>::write(const Parameters_ & params) const {
  Log::trace() << "GeoVaLs<OBS>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  gvals_->write(params);
  Log::trace() << "GeoVaLs<OBS>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void GeoVaLs<OBS>::print(std::ostream & os) const {
  Log::trace() << "GeoVaLs<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *gvals_;
  Log::trace() << "GeoVaLs<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_GEOVALS_H_
