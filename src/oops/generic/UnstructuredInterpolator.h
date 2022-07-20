/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Metadata.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL>
class UnstructuredInterpolator : public util::Printable,
                                 private util::ObjectCounter<UnstructuredInterpolator<MODEL>> {
  typedef State<MODEL>     State_;
  typedef Geometry<MODEL>  Geometry_;
  typedef Increment<MODEL> Increment_;

 public:
  static const std::string classname() {return "oops::UnstructuredInterpolator";}

  UnstructuredInterpolator(const eckit::Configuration &, const Geometry_ &,
                           const std::vector<double> &, const std::vector<double> &);

  void apply(const Variables &, const State_ &, const std::vector<bool> &,
             std::vector<double> &) const;
  void apply(const Variables &, const Increment_ &, const std::vector<bool> &,
             std::vector<double> &) const;
  void applyAD(const Variables &, Increment_ &, const std::vector<bool> &,
               const std::vector<double> &) const;

 private:
  void apply(const Variables &, const atlas::FieldSet &, const std::vector<bool> &,
             std::vector<double> &) const;

  void applyPerLevel(const std::string &, const std::vector<bool> &,
                     const atlas::array::ArrayView<double, 2> &,
                     std::vector<double>::iterator &, const size_t &) const;
  void applyPerLevelAD(const std::string &, const std::vector<bool> &,
                       atlas::array::ArrayView<double, 2> &,
                       std::vector<double>::const_iterator &, const size_t &) const;
  void print(std::ostream &) const override;

  std::string interp_method_;
  int nninterp_;
  size_t nout_;
  std::vector<std::vector<size_t>> interp_i_;
  std::vector<std::vector<double>> interp_w_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
UnstructuredInterpolator<MODEL>::UnstructuredInterpolator(const eckit::Configuration & config,
                                                          const Geometry_ & grid,
                                                          const std::vector<double> & lats_out,
                                                          const std::vector<double> & lons_out)
  : interp_method_(), nninterp_(0), nout_(0), interp_i_(), interp_w_()
{
  Log::trace() << "UnstructuredInterpolator::UnstructuredInterpolator start" << std::endl;
  util::Timer timer("oops::UnstructuredInterpolator", "UnstructuredInterpolator");

  const atlas::Geometry earth(atlas::util::Earth::radius());
  const double close = 1.0e-10;

  // This is a new option for this class, so isn't in any YAMLs yet!
  interp_method_ = config.getString("interpolation method", "barycentric");
  ASSERT(interp_method_ == "barycentric" || interp_method_ == "inverse distance");

  // Compute weights
  ASSERT(lats_out.size() == lons_out.size());
  nout_ = lats_out.size();
  nninterp_ = config.getInt("nnearest", 4);
  interp_i_.resize(nout_, std::vector<size_t>(nninterp_));
  interp_w_.resize(nout_, std::vector<double>(nninterp_, 0.0));

  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    const atlas::util::KDTree<size_t>::ValueList neighbours =
                          grid.closestPoints(lats_out[jloc], lons_out[jloc], nninterp_);

    // Barycentric and inverse-distance interpolation both rely on indices, 1/distances
    size_t jj = 0;
    for (const atlas::util::KDTree<size_t>::Value & val : neighbours) {
      interp_i_[jloc][jj] = val.payload();
      interp_w_[jloc][jj] = 1.0 / val.distance();
      ++jj;
    }
    ASSERT(jj == nninterp_);

    if (interp_w_[jloc][0] > 1e10) {
      // Handle edge case where output point is close to one input point => interp_w very large
      // Atlas returns the neighbors in nearest-first order, so only need to check first element
      for (size_t jn = 0; jn < nninterp_; ++jn) interp_w_[jloc][jn] = 0.0;
      interp_w_[jloc][0] = 1.0;
    } else if (interp_method_ == "barycentric") {
      // Barycentric weights
      std::vector<double> bw(nninterp_);
      for (size_t j = 0; j < nninterp_; ++j) {
        double wprod = 1.0;
        const auto& p1 = neighbours[j].point();
        for (size_t k = 0; k < nninterp_; ++k) {
          if (k != j) {
            const auto& p2 = neighbours[k].point();
            const double dist = earth.distance(p1, p2);
            wprod *= std::max(dist, close);
          }
        }
        bw[j] = 1.0 / wprod;
      }
      // Interpolation weights from barycentric weights
      double wsum = 0.0;
      for (size_t j = 0; j < nninterp_; ++j) wsum += interp_w_[jloc][j] * bw[j];
      for (size_t j = 0; j < nninterp_; ++j) {
        interp_w_[jloc][j] *= bw[j] / wsum;
        ASSERT(interp_w_[jloc][j] >= 0.0 && interp_w_[jloc][j] <= 1.0);
      }
    } else {
      // Inverse-distance interpolation weights
      double wsum = 0.0;
      for (size_t j = 0; j < nninterp_; ++j) wsum += interp_w_[jloc][j];
      for (size_t j = 0; j < nninterp_; ++j) {
        interp_w_[jloc][j] /= wsum;
        ASSERT(interp_w_[jloc][j] >= 0.0 && interp_w_[jloc][j] <= 1.0);
      }
    }
  }
  Log::trace() << "UnstructuredInterpolator::UnstructuredInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator<MODEL>::apply(const Variables & vars, const State_ & xx,
                                            const std::vector<bool> & mask,
                                            std::vector<double> & locvals) const
{
  this->apply(vars, xx.fieldSet(), mask, locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator<MODEL>::apply(const Variables & vars, const Increment_ & dx,
                                            const std::vector<bool> & mask,
                                            std::vector<double> & locvals) const
{
  this->apply(vars, dx.fieldSet(), mask, locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator<MODEL>::apply(const Variables & vars, const atlas::FieldSet & fset,
                                            const std::vector<bool> & mask,
                                            std::vector<double> & vals) const {
  Log::trace() << "UnstructuredInterpolator::apply starting" << std::endl;
  util::Timer timer("oops::UnstructuredInterpolator", "apply");

  size_t nflds = 0;
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf];
    nflds += fset.field(fname).levels();
  }
  vals.resize(nout_ * nflds);

  auto current = vals.begin();
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf];
    atlas::Field fld = fset.field(fname);

    const std::string interp_type = fld.metadata().get<std::string>("interp_type");
    ASSERT(interp_type == "default" || interp_type == "integer" || interp_type == "nearest");

    const atlas::array::ArrayView<double, 2> fldin = atlas::array::make_view<double, 2>(fld);
    for (size_t jlev = 0; jlev < fldin.shape(1); ++jlev) {
      this->applyPerLevel(interp_type, mask, fldin, current, jlev);
    }
  }
  Log::trace() << "UnstructuredInterpolator::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator<MODEL>::applyAD(const Variables & vars, Increment_ & dx,
                                              const std::vector<bool> & mask,
                                              const std::vector<double> & vals) const {
  Log::trace() << "UnstructuredInterpolator::applyAD starting" << std::endl;
  util::Timer timer("oops::UnstructuredInterpolator", "applyAD");

  std::vector<double>::const_iterator current = vals.begin();
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf];
    atlas::Field fld = dx.fieldSet().field(fname);

//    const std::string interp_type = fld.metadata().get<std::string>("interp_type");
//    ASSERT(interp_type == "default" || interp_type == "integer" || interp_type == "nearest");
    const std::string interp_type = "default";

    atlas::array::ArrayView<double, 2> fldin = atlas::array::make_view<double, 2>(fld);
    for (size_t jlev = 0; jlev < fldin.shape(1); ++jlev) {
      this->applyPerLevelAD(interp_type, mask, fldin, current, jlev);
    }
  }
  Log::trace() << "UnstructuredInterpolator::applyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator<MODEL>::applyPerLevel(
    const std::string & interp_type, const std::vector<bool> & mask,
    const atlas::array::ArrayView<double, 2> & gridin,
    std::vector<double>::iterator & gridout, const size_t & ilev) const {
  ASSERT(mask.size() == nout_);
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    if (mask[jloc]) {
      *gridout = 0.0;
      if (interp_type == "default") {
        for (size_t jj = 0; jj < nninterp_; ++jj) {
          *gridout += interp_w_[jloc][jj] * gridin(interp_i_[jloc][jj], ilev);
        }
      } else if (interp_type == "integer") {
        // Find which integer value has largest weight in the stencil. We do this by taking two
        // passes through the (usually short) data: first to identify range of values, then to
        // determine weights for each integer.
        // Note that a std::map would be shorter to code, because it would avoid needing to find
        // the range of possible integer values, but vectors are almost always much more efficient.
        int minval = std::numeric_limits<int>().max();
        int maxval = std::numeric_limits<int>().min();
        for (size_t jj = 0; jj < nninterp_; ++jj) {
         minval = std::min(minval, static_cast<int>(std::round(gridin(interp_i_[jloc][jj], ilev))));
         maxval = std::max(maxval, static_cast<int>(std::round(gridin(interp_i_[jloc][jj], ilev))));
        }
        std::vector<double> int_weights(maxval - minval + 1, 0.0);
        for (size_t jj = 0; jj < nninterp_; ++jj) {
          const int this_int = std::round(gridin(interp_i_[jloc][jj], ilev));
          int_weights[this_int - minval] += interp_w_[jloc][jj];
        }
        *gridout = minval + std::distance(int_weights.begin(),
            std::max_element(int_weights.begin(), int_weights.end()));
      } else if (interp_type == "nearest") {
        *gridout = gridin(interp_i_[jloc][0], ilev);
      } else {
        throw eckit::BadValue("Unknown interpolation type");
      }
    }
    ++gridout;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator<MODEL>::applyPerLevelAD(const std::string & interp_type,
                                                      const std::vector<bool> & mask,
                                                      atlas::array::ArrayView<double, 2> & gridin,
                                                      std::vector<double>::const_iterator & gridout,
                                                      const size_t & ilev) const {
  ASSERT(mask.size() == nout_);
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    if (mask[jloc]) {
      if (interp_type == "default") {
        for (size_t jj = 0; jj < nninterp_; ++jj) {
          gridin(interp_i_[jloc][jj], ilev) += interp_w_[jloc][jj] * *gridout;
        }
      } else if (interp_type == "integer") {
        throw eckit::BadValue("No adjoint for integer interpolation");
      } else if (interp_type == "nearest") {
        gridin(interp_i_[jloc][0], ilev) += *gridout;
      } else {
        throw eckit::BadValue("Unknown interpolation type");
      }
    }
    ++gridout;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator<MODEL>::print(std::ostream & os) const
{
  os << "UnstructuredInterpolator<" << MODEL::name() << ">";
}

// -----------------------------------------------------------------------------

}  // namespace oops
