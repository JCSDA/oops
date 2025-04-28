/*
 * (C) Crown Copyright 2025, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/assimilation/SubensembleSplitter.h"

#include <Eigen/Dense>
#include <algorithm>
#include <random>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

SubensembleSplitter::SubensembleSplitter(const size_t nens,
                                         const eckit::LocalConfiguration & config)
    : nens_(nens), neig_(1)
{
    ctorHelper(config);
}

// -----------------------------------------------------------------------------

SubensembleSplitter::SubensembleSplitter(const size_t nens, const size_t neig,
                                         const eckit::LocalConfiguration & config)
    : nens_(nens), neig_(neig)
{
    ctorHelper(config);
}

// -----------------------------------------------------------------------------
void SubensembleSplitter::ctorHelper(const eckit::LocalConfiguration & config) {
    // Reading from config
    nsubens_ = config.getUnsigned("number of subensembles");
    splittingMethod_ = config.getString("splitting method", "Contiguous");
    seed_ = config.getUnsigned("random seed", 0);

    // Validation
    if ((nsubens_ < 1) || (nsubens_ > nens_)) {
      oops::Log::error() << "nsubens = " << nsubens_
                         << " out of bounds, must be in the closed interval [1, "
                         << nens_ << "]\n";
      throw eckit::BadParameter("nsubens out of bounds, must be in the closed interval [1, nens]");
    }
    if (nens_ % nsubens_ != 0) {
      oops::Log::error() << "nens = " << nens_
                         << " not divisible by nsubens = "
                         << nsubens_ << "\n";
      throw eckit::BadParameter("nens not divisible by nsubens");
    }
    if (neig_ < 1) {
      oops::Log::error() << "neig cannot be less than 1"
                         << " but was found to be " << neig_ << "\n";
      throw eckit::BadParameter("neig cannot be less than 1");
    }

    // Private variable initialisation
    subens_.resize(nens_);
    subsize_ = nens_/nsubens_;
}

// -----------------------------------------------------------------------------

void SubensembleSplitter::split() {
  std::vector<size_t> indices(nens_);
  std::iota(indices.begin(), indices.end(), 0);
  if (splittingMethod_ == "Contiguous") {
    contiguousPartition(indices);
  } else if (splittingMethod_ == "Random") {
    std::default_random_engine rng;
    rng.seed(seed_);
    std::shuffle(std::begin(indices), std::end(indices), rng);
    contiguousPartition(indices);
  } else {
    oops::Log::error() << "Splitting method '" << splittingMethod_ << "' not implemented\n";
    throw eckit::BadParameter("Splitting method not implemented");
  }
}

// -----------------------------------------------------------------------------

void SubensembleSplitter::contiguousPartition(const std::vector<size_t> & indices) {
  for (size_t subens = 0; subens < nsubens_; ++subens) {
    const size_t loop_start = subens*subsize_;
    const size_t loop_end = (subens+1)*subsize_;
    for (size_t part = loop_start; part < loop_end; ++part) {
      const size_t index = indices[part];
      subens_[index] = subens;
    }
  }
}

// -----------------------------------------------------------------------------

const std::tuple<Eigen::MatrixXf, std::vector<size_t>>
      SubensembleSplitter::getExclusionTensors(const size_t xclsub, const bool modulated) {
  // Validating excluded subensemble
  if ((xclsub < 0) || (xclsub >= nsubens_)) {
    oops::Log::error() << "Excluded subensemble " << xclsub
                       << " is not in the valid interval [0, "
                       << nsubens_ << ")\n";
    throw eckit::BadValue("Excluded subensemble is not in "
                          "the valid interval [0, nsubens)");
  }

  // Setting up matrix sizes
  size_t neig;
  if (modulated) {
    neig = neig_;
  } else {
    neig = 1;
  }
  const size_t nmat = nens_*neig;
  const size_t nhat = nmat - subsize_*neig;
  Eigen::MatrixXf proj = Eigen::MatrixXf::Zero(nmat, nhat);
  std::vector<size_t> excluded;

  // Building projection matrix
  size_t col = 0;
  for (size_t iens = 0; iens < nens_; ++iens) {
    const bool isExcluded = (subens_[iens] == xclsub);
    if (!isExcluded) {
      for (size_t ieig = 0; ieig < neig; ++ieig) {
          const size_t row = iens*neig + ieig;
          proj(row, col) = 1.0f;
          ++col;
      }
    } else {
        excluded.emplace_back(iens);
    }
  }
  excluded.shrink_to_fit();
  return std::make_tuple(proj, excluded);
}

// -----------------------------------------------------------------------------

}  // namespace oops
