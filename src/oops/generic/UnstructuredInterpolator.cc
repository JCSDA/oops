/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/UnstructuredInterpolator.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Point.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Timer.h"

namespace detail {
template <size_t N>
std::array<size_t, N> decreasing_permutation(const std::array<double, N> & arr) {
  std::array<std::size_t, N> p{};
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(), [&arr](const size_t i, const size_t j){ return arr[i] > arr[j]; });
  return p;
}

template <typename T, size_t N>
std::array<T, N> permute_array(const std::array<T, N> & arr, const std::array<std::size_t, N> & p) {
  std::array<T, N> sorted{};
  std::transform(p.begin(), p.end(), sorted.begin(), [&](const size_t i){ return arr[i]; });
  return sorted;
}
}  // namespace detail

namespace oops {

// -----------------------------------------------------------------------------

UnstructuredInterpolator::UnstructuredInterpolator(const eckit::Configuration & config,
                                                   const GeometryData & geomData,
                                                   const std::vector<double> & lats_out,
                                                   const std::vector<double> & lons_out)
  : atlasbase::Interpolator(config, geomData, lats_out, lons_out),
    geom_(geomData), nout_(0), interp_matrices_{}
{
  Log::trace() << "UnstructuredInterpolator::UnstructuredInterpolator start" << std::endl;
  util::Timer timer("oops::UnstructuredInterpolator", "UnstructuredInterpolator");

  ASSERT(lats_out.size() == lons_out.size());
  nout_ = lats_out.size();

  regionalNnFillDistance_ = config.getDouble("regional nn fill distance in km", 0.0) * 1000.0;
  enableRegionalCheck_ = config.getBool("regional check enabled", true);

  computeUnmaskedInterpMatrix(lats_out, lons_out);

  Log::trace() << "UnstructuredInterpolator::UnstructuredInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::apply(const Variables & vars, const atlas::FieldSet & fset,
                                     const std::vector<bool> & target_mask,
                                     std::vector<double> & vals) const {
  Log::trace() << "UnstructuredInterpolator::apply starting" << std::endl;
  util::Timer timer("oops::UnstructuredInterpolator", "apply");

  if (nout_ == 0) { return; }
  ASSERT(target_mask.size() == nout_);

  size_t nflds = 0;
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf].name();
    nflds += fset.field(fname).shape(1);
  }
  vals.resize(nout_ * nflds);

  auto current = vals.begin();
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf].name();
    atlas::Field & fld = fset.field(fname);  // const in principle, but intel can't compile that

    const std::string interp_type = fld.metadata().get<std::string>("interp_type");
    ASSERT(interp_type == "default" || interp_type == "integer" || interp_type == "nearest");

    // Mask is optional -- no metadata signals unmasked interpolation
    // Warning: if the model code typoes the name of the metadata field "mask",
    // then the code below will silently skip the masking and proceed with unmasked interpolation.
    // Requiring the mask metadata to be always present would increase robustness, but would require
    // all models to adapt.
    std::string maskName = unmaskedName_;
    if (fld.metadata().has("mask")) {
      maskName = fld.metadata().get<std::string>("mask");
      ASSERT(geom_.has(maskName));
      const atlas::Field & source_mask_fld = geom_.getField(maskName);
      ASSERT(source_mask_fld.shape(0) == fld.shape(0));
      ASSERT(source_mask_fld.shape(1) == 1);  // For now, support 2D masks only
      const auto source_mask = atlas::array::make_view<double, 2>(source_mask_fld);

      // Compute the masked interpolation matrix for this mask, if not previously done
      if (interp_matrices_.find(maskName) == interp_matrices_.end()) {
        computeMaskedInterpMatrix(maskName, source_mask);
      }
    }

    // Get interpolation matrix for this mask
    const auto & interpMatrix = interp_matrices_.at(maskName);

    // Get array views into source Field and target vector
    const auto source = atlas::array::make_view<double, 2>(fld);
    const auto vals_shape =
        atlas::array::ArrayShape{static_cast<int>(nout_), source.shape(1)};
    std::unique_ptr<atlas::array::Array> vals_array(
        atlas::array::Array::wrap<double>(&*current, vals_shape));
    auto target = atlas::array::make_view<double, 2>(*vals_array);

    // Interpolate
    this->doApply(interpMatrix, interp_type, target_mask, source, target);

    // Advance iterator to next variable without vals
    const int vals_step = nout_ * source.shape(1);
    std::advance(current, vals_step);
  }
  Log::trace() << "UnstructuredInterpolator::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::applyAD(const Variables & vars, atlas::FieldSet & fset,
                                       const std::vector<bool> & target_mask,
                                       const std::vector<double> & vals) const {
  Log::trace() << "UnstructuredInterpolator::applyAD starting" << std::endl;
  util::Timer timer("oops::UnstructuredInterpolator", "applyAD");

  if (nout_ == 0) { return; }
  ASSERT(target_mask.size() == nout_);

  std::vector<double>::const_iterator current = vals.begin();
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf].name();
    atlas::Field & fld = fset.field(fname);

//    const std::string interp_type = fld.metadata().get<std::string>("interp_type");
//    ASSERT(interp_type == "default" || interp_type == "integer" || interp_type == "nearest");
    const std::string interp_type = "default";

    // Mask is optional -- no metadata signals unmasked interpolation
    std::string maskName = unmaskedName_;
    if (fld.metadata().has("mask")) {
      maskName = fld.metadata().get<std::string>("mask");
      ASSERT(geom_.has(maskName));
      const atlas::Field & source_mask_fld = geom_.getField(maskName);
      ASSERT(source_mask_fld.shape(0) == fld.shape(0));
      ASSERT(source_mask_fld.shape(1) == 1);  // For now, support 2D masks only
      const auto source_mask = atlas::array::make_view<double, 2>(source_mask_fld);

      // Compute the masked interpolation matrix for this mask, if not previously done
      if (interp_matrices_.find(maskName) == interp_matrices_.end()) {
        computeMaskedInterpMatrix(maskName, source_mask);
      }
    }

    // Get interpolation matrix for this mask
    const auto & interpMatrix = interp_matrices_.at(maskName);

    // Get array views into source Field and target vector
    auto source = atlas::array::make_view<double, 2>(fld);
    const auto vals_shape =
        atlas::array::ArrayShape{static_cast<int>(nout_), source.shape(1)};
    // Horrible cast to allow a standard ArrayView<double> from a const std::vector
    const std::unique_ptr<atlas::array::Array> vals_array(
        atlas::array::Array::wrap<double>(const_cast<double*>(&*current), vals_shape));
    const auto target = atlas::array::make_view<double, 2>(*vals_array);

    // Interpolate
    this->doApplyAD(interpMatrix, interp_type, target_mask, source, target);

    // Advance iterator to next variable without vals
    const int vals_step = nout_ * source.shape(1);
    std::advance(current, vals_step);
  }
  Log::trace() << "UnstructuredInterpolator::applyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::doApply(
    const InterpMatrix & interpMatrix,
    const std::string & interp_type,
    const std::vector<bool> & target_mask,
    const atlas::array::ArrayView<double, 2> & source,
    atlas::array::ArrayView<double, 2> & target) const {
  const int nb_levels = source.shape(1);
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    if (target_mask[jloc]) {
      // Set the location column to zero before interpolating.
      // If there is no valid stencil for this location, set to missing.
      if (interpMatrix.targetHasValidStencil[jloc]) {
        for (int jlev = 0; jlev < nb_levels; ++jlev) {
          target(jloc, jlev) = 0.;
        }
      } else {
        for (int jlev = 0; jlev < nb_levels; ++jlev) {
          target(jloc, jlev) = util::missingValue<double>();
        }
        continue;  // invalid stencil, so exit early
      }

      const std::vector<size_t> & interp_is = interpMatrix.stencils[jloc];
      const std::vector<double> & interp_ws = interpMatrix.weights[jloc];

      if (interp_type == "default") {
        for (size_t jj = 0; jj < nstencil_; ++jj) {
          for (int jlev = 0; jlev < nb_levels; ++jlev) {
            target(jloc, jlev) += interp_ws[jj] * source(interp_is[jj], jlev);
          }
        }
      } else if (interp_type == "nearest") {
        // Return value from closest unmasked source point
        for (size_t jj = 0; jj < nstencil_; ++jj) {
          if (interp_ws[jj] > 1.0e-9) {  // use a small tolerance to allow for roundoff in weights
            for (int jlev = 0; jlev < nb_levels; ++jlev) {
              target(jloc, jlev) = source(interp_is[jj], jlev);
            }
            break;
          }
        }
      } else if (interp_type == "integer") {
        // Use with caution: this is a slow and nonlinear "voting" scheme.
        // The scheme finds which integer value has largest weight across the interpolation
        // stencil. This is done by taking two passes through the (usually short) data: first to
        // identify the range of the exisitng values, then to determine weights for each integer.
        // TODO(core team): this is fv3-jedi specific code that needs to be removed from oops.
        for (int jlev = 0; jlev < nb_levels; ++jlev) {
          int minval = std::numeric_limits<int>().max();
          int maxval = std::numeric_limits<int>().min();
          for (size_t jj = 0; jj < nstencil_; ++jj) {
            minval = std::min(minval, static_cast<int>(std::round(source(interp_is[jj], jlev))));
            maxval = std::max(maxval, static_cast<int>(std::round(source(interp_is[jj], jlev))));
          }
          std::vector<double> int_weights(maxval - minval + 1, 0.0);
          for (size_t jj = 0; jj < nstencil_; ++jj) {
            const int this_int = std::round(source(interp_is[jj], jlev));
            int_weights[this_int - minval] += interp_ws[jj];
          }
          target(jloc, jlev) = minval + std::distance(int_weights.begin(),
              std::max_element(int_weights.begin(), int_weights.end()));
        }
      } else {
        throw eckit::BadValue("Unknown interpolation type");
      }
    }
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::doApplyAD(
    const InterpMatrix & interpMatrix,
    const std::string & interp_type,
    const std::vector<bool> & target_mask,
    atlas::array::ArrayView<double, 2> & source,
    const atlas::array::ArrayView<double, 2> & target) const {
  const int nb_levels = source.shape(1);
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    if (target_mask[jloc]) {
      // (Adjoint of) If there is no valid stencil for this location, set to missing.
      if (!interpMatrix.targetHasValidStencil[jloc]) {
        continue;
      }

      const std::vector<size_t> & interp_is = interpMatrix.stencils[jloc];
      const std::vector<double> & interp_ws = interpMatrix.weights[jloc];

      if (interp_type == "default") {
        for (size_t jj = 0; jj < nstencil_; ++jj) {
          for (int jlev = 0; jlev < nb_levels; ++jlev) {
            source(interp_is[jj], jlev) += interp_ws[jj] * target(jloc, jlev);
          }
        }
      } else if (interp_type == "nearest") {
        // (Adjoint of) Return value from closest unmasked source point
        for (size_t jj = 0; jj < nstencil_; ++jj) {
          if (interp_ws[jj] > 1.0e-9) {  // use a small tolerance to allow for roundoff in weights
            for (int jlev = 0; jlev < nb_levels; ++jlev) {
              source(interp_is[jj], jlev) += target(jloc, jlev);
            }
            break;
          }
        }
      } else if (interp_type == "integer") {
        throw eckit::BadValue("No adjoint for integer interpolation");
      } else {
        throw eckit::BadValue("Unknown interpolation type");
      }
    }
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::computeUnmaskedInterpMatrix(
    std::vector<double> lats_out,
    std::vector<double> lons_out) const {
  // Check matrix hasn't already been set
  ASSERT(interp_matrices_.find(unmaskedName_) == interp_matrices_.end());

  // Compute interpolation matrix with no source-point mask
  interp_matrices_.insert(
      std::make_pair(unmaskedName_, InterpMatrix{
        std::vector<bool>(nout_, true),
        std::vector<std::vector<size_t>>(nout_, std::vector<size_t>(nstencil_)),
        std::vector<std::vector<double>>(nout_, std::vector<double>(nstencil_, 0.0))}));

  // The logic switch on enableRegionalCheck_ enables an "unsafe" mode where the
  // atlas-based isRegional() call is skipped, thereby not ensuring the regional
  // nn fill feature is only activated in a regional domain. This is unsafe
  // because using this feature in a global domain could mask real geometry
  // construction errors... so ideally it would be disabled.
  // The rationale for providing the unsafe mode is to temporarily support
  // models whose atlas geometry setup backing isRegional() may be buggy.
  // TODO(FH): Remove enableRegionalCheck_ and unsafe mode, by instead ensuring
  //           all models are able to provide an accurate isRegional()
  bool enableRegionalNnFill = false;
  if (enableRegionalCheck_) {
    const bool isRegional = util::isRegional(geom_.functionSpace(), geom_.fieldSet());
    enableRegionalNnFill = (isRegional && regionalNnFillDistance_ > 0.0);
  } else {
    enableRegionalNnFill = true;
  }

  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    std::array<int, 3> indices{};
    std::array<double, 3> baryCoords{};
    const bool validTriangle = geom_.containingTriangleAndBarycentricCoords(
        lats_out[jloc], lons_out[jloc], indices, baryCoords);

    // Edge case: target point outside of source grid, can occur for regional models
    if (!validTriangle) {
      if (enableRegionalNnFill) {
        const auto search_result = geom_.closestPointWithinRadius(lats_out[jloc], lons_out[jloc],
                                                                  regionalNnFillDistance_);
        if (search_result.has_value()) {
          // NN extrapolation: full weight on the single nearest source point
          const int index = search_result.value();
          auto & interp_is = interp_matrices_.at(unmaskedName_).stencils[jloc];
          auto & interp_ws = interp_matrices_.at(unmaskedName_).weights[jloc];
          interp_is[0] = static_cast<size_t>(index);
          interp_is[1] = static_cast<size_t>(index);
          interp_is[2] = static_cast<size_t>(index);
          interp_ws[0] = 1.0;
          interp_ws[1] = 0.0;
          interp_ws[2] = 0.0;
          // targetHasValidStencil[jloc] remains true (as initialized)
          continue;
        }
      }
      interp_matrices_[unmaskedName_].targetHasValidStencil[jloc] = false;
      continue;
    }

    // Reorder points from nearest to furthest (barycentric coords from largest to smallest)
    // This ordering of the coefficients is used in nearest-neighbor interpolations
    const auto permutation = detail::decreasing_permutation(baryCoords);
    indices = detail::permute_array(indices, permutation);
    baryCoords = detail::permute_array(baryCoords, permutation);

    double wsum = 0.0;
    for (size_t j = 0; j < nstencil_; ++j) {
      wsum += baryCoords[j];
      ASSERT(baryCoords[j] >= 0.0 && baryCoords[j] <= 1.0);
    }
    ASSERT(fabs(wsum - 1.0) < 1e-14);

    // Store indices and weights into InterpMatrix datastructure
    std::vector<size_t> & interp_is = interp_matrices_.at(unmaskedName_).stencils[jloc];
    std::vector<double> & interp_ws = interp_matrices_.at(unmaskedName_).weights[jloc];
    for (size_t j = 0; j < nstencil_; ++j) {
      interp_is[j] = indices[j];
      interp_ws[j] = baryCoords[j];
    }
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::computeMaskedInterpMatrix(
    const std::string & maskName,
    const atlas::array::ArrayView<double, 2> & source_mask) const
{
  // Check unmasked matrix is already computed
  ASSERT(interp_matrices_.find(unmaskedName_) != interp_matrices_.end());
  // Check matrix hasn't already been computed for this mask
  ASSERT(interp_matrices_.find(maskName) == interp_matrices_.end());

  // Copy unmasked matrix, then modify it below
  interp_matrices_[maskName] = interp_matrices_[unmaskedName_];

  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    // Edge case: unmasked interp stencil is already invalid => masked matrix also invalid
    if (!interp_matrices_[unmaskedName_].targetHasValidStencil[jloc]) {
      interp_matrices_[maskName].targetHasValidStencil[jloc] = false;
      continue;
    }
    std::vector<size_t> & interp_is = interp_matrices_[maskName].stencils[jloc];
    std::vector<double> & interp_ws = interp_matrices_[maskName].weights[jloc];

    // Sum up mask weights, will be used to renormalize interpolation weights
    double normalization = 0.0;
    for (size_t jj = 0; jj < nstencil_; ++jj) {
      ASSERT(source_mask(interp_is[jj], 0) >= 0.0 && source_mask(interp_is[jj], 0) <= 1.0);
      normalization += interp_ws[jj] * source_mask(interp_is[jj], 0);
    }

    if (normalization <= 1e-9) {
      // Edge case: all source points are masked out, so can't interpolate to this target point
      interp_matrices_[maskName].targetHasValidStencil[jloc] = false;
    } else {
      // Standard case: renormalize
      for (size_t jj = 0; jj < nstencil_; ++jj) {
        interp_ws[jj] *= source_mask(interp_is[jj], 0) / normalization;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::print(std::ostream & os) const { os << classname(); }

// -----------------------------------------------------------------------------

}  // namespace oops
