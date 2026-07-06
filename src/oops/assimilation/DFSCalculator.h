/*
 * (C) Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_DFSCALCULATOR_H_
#define OOPS_ASSIMILATION_DFSCALCULATOR_H_

#include <Eigen/Dense>

#include <cstddef>
#include <iomanip>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Departures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"

namespace oops {

/**
 * @brief Degrees of Freedom for Signal (DFS) Calculator for LETKF.
 *
 * This class calculates the DFS per observation and save then into obspace. 
 * It also writes the average DFS per type to the log file.
 *
 * Workflow:
 *   1) For each grid point, call accumulate(...) with:
 *      - omb: defines the global serial layout and block structure.
 *      - locvector: selects local/assimilated obs.
 *      - Yb: ensemble perturbations in obs space (nens x nloc).
 *      - Wa: transform matrix that maps Yb→Ya (nens x nens).
 *      - invR: inverse obs error variance for local obs (length nloc).
 *   2) After looping over all grid points, call finalize() to average
 *      contributions for observations used multiple times.
 *   3) Call computeBlockStats(...) to compute the DFS per type and 
 *      call printBlockStats(...) to save the result in the log file.
 *   4) Write DFS per observation in obspace via savePerObservation(...).
 *
 * Notes: DFS per local obs is computed as diag( (Wa*Yb)ᵀ (Wa*Yb) ) ⊙ invR / (nens-1).
 * 
 * Reference:
 * Hu, G., Dance, S. L., Fowler, A., Simonin, D., & Waller, J. (2025). 
 * Assessing the influence of observations in ensemble-based data assimilation systems. 
 * Journal of Advances in Modeling Earth Systems, 17, e2024MS004809. 
 * https://doi.org/10.1029/2024MS004809
 */

template <typename OBS>
class DFSCalculator {
 public:
  using index_t     = Eigen::Index;
  using Departures_ = oops::Departures<OBS>;
  using ObsSpaces_  = oops::ObsSpaces<OBS>;

  /// Construct accumulators sized to total number of observations.
  explicit DFSCalculator(std::size_t nobs)
      : dfs_(Eigen::VectorXd::Zero(static_cast<index_t>(nobs))),
        usage_(Eigen::VectorXd::Zero(static_cast<index_t>(nobs))) {}
  /// Simple per-observation block summary (label, count, mean).
  struct BlockStat {
    std::string   label;      // e.g., "Stream", "Wind", "WSpeed"
    std::size_t   nobs{0};    // selected obs count in this block
    double        mean{std::numeric_limits<double>::quiet_NaN()};
  };

  /// Compute local DFS contributions for a single grid point and
  /// accumulate them into the global observation vector.
  void accumulate(const Departures_ & omb,                                  // defines global order
                  const Departures_ & locvector,                            // local selection
                  const Eigen::Ref<const Eigen::MatrixXf> & Yb,             // (nens x nloc)
                  const Eigen::Ref<const Eigen::MatrixXf> & Wa,             // (nens x nens)
                  const Eigen::Ref<const Eigen::VectorXd> & invR) {         // (nloc)
    const index_t nloc = Yb.cols();
    const index_t nens = Yb.rows();

    if (nloc == 0) return;  // skip the grid point if no observations are used
    ASSERT(nens > 1);  // avoid division by zero in (nens-1)
    const double invNe1 = 1.0 / static_cast<double>(nens - 1);
    ensureBuffers(Wa.rows(), nloc);  // ensure scratch buffers are the right shape.

    // Analysis ensemble perturbations in observation space
    Ya_buf_.noalias() = Wa * Yb;

    // Local DFS values in the same order as packEigen(locvector)
    diagS_buf_.noalias() =
      Ya_buf_.cast<double>().colwise().squaredNorm().transpose().cwiseProduct(invR);
    diagS_buf_ *= invNe1;

    // Map local indices → global indices.
    // gid_blocks[b] contains serial indices *within block b* for selected obs.
    const auto gid_blocks = omb.maskAndSerialIndices(locvector);
    // block_off[b] is the global offset where block b starts in the concatenated layout.
    const auto &block_off = blockOffsetsCached(omb);
    // Accumulate local DFS contributions into the global vector.
    index_t j = 0;
    for (std::size_t b = 0; b < gid_blocks.size(); ++b) {
      const std::size_t off = block_off[b];
      for (std::size_t g_local : gid_blocks[b]) {
        const std::size_t g = off + g_local;  // convert local id → global id
        const index_t gid = static_cast<index_t>(g);
        dfs_(gid)   += diagS_buf_(j);  // sum contributions across grid points
        usage_(gid) += 1.0;            // count how many times this obs contributed
        ++j;
      }
    }
  }

  /// Average contributions; call once when done accumulating
  void finalize() {
    if (dfs_.size() == 0) return;
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    dfs_.array() = (usage_.array() > 0.0).select(dfs_.array() / usage_.array(), NaN);
  }

  /// Compute average DFS per observation block.
  std::vector<BlockStat>
  computeBlockStats(const Departures_ & omb,
                    const Departures_ & mask,
                    const eckit::Configuration & observersConf) const {
    const auto gid_blocks = omb.maskAndSerialIndices(mask);
    auto labels = extractObsTypes(observersConf);
    const auto &block_off = blockOffsetsCached(omb);

    std::vector<BlockStat> stats;
    stats.reserve(gid_blocks.size());

    for (std::size_t b = 0; b < gid_blocks.size(); ++b) {
      const auto & local_ids = gid_blocks[b];
      const auto  off  = block_off[b];
      const auto  nb   = local_ids.size();

      BlockStat s;
      s.label = (b < labels.size() ? labels[b] : std::string());

      double sum = 0.0;
      if (nb > 0) {
        for (const std::size_t gl : local_ids) {
          const index_t gid = static_cast<index_t>(off + gl);
          sum += dfs_(gid);
        }
      }

      int globalNb = static_cast<int>(nb);

      s.nobs = static_cast<std::size_t>(globalNb);
      if (globalNb > 0) {
        s.mean = sum / static_cast<double>(globalNb);
      }
      stats.push_back(std::move(s));
    }
    return stats;
  }

  /// Print DFS per block to the oops logger.
  void printBlockStats(const std::vector<BlockStat> & stats) const {
    for (const auto & s : stats) {
      oops::Log::info()
        << "DFS "
        << std::left  << std::setw(24) << s.label
        << " nobs="   << std::right << std::setw(6) << s.nobs
        << " mean="   << std::scientific << std::setprecision(6) << s.mean
        << std::defaultfloat
        << std::endl;
    }
  }
  /// Scatter finalized DFS into each ObsSpace and save into the nc file.
  void savePerObservation(const Departures_ &omb,
                          const Departures_ &mask,
                          const ObsSpaces_  &obspaces,
                          const std::string &varName) const {
    // packed indices (same order as packEigen(mask)) and cached block offsets
    const auto gid_blocks = omb.maskAndSerialIndices(mask);
    const auto &block_off = blockOffsetsCached(omb);

    // container that mirrors ObsSpaces for output
    Departures_ dfsDep(obspaces);
    const Eigen::VectorXd &dfsGlobal = this->getDFS();

    for (std::size_t b = 0; b < gid_blocks.size(); ++b) {
      const std::size_t off = block_off[b];

      // Serialize full raw length into a plain vector<double>
      std::vector<double> raw;
      dfsDep[b].serialize(raw);

      // inds maps packed positions → raw indices within this block
      const auto &inds = dfsDep[b].maskAndSerialIndices(mask[b]);
      for (std::size_t k = 0; k < inds.size(); ++k) {
        const index_t g = static_cast<index_t>(off + gid_blocks[b][k]);
        raw[inds[k]] = dfsGlobal(g);
      }

      // Push the updated raw data back into dfsDep for this block
      std::size_t ind = 0;
      dfsDep[b].deserialize(raw, ind);
    }
    dfsDep.save(varName);
  }

  /// Accessors
  const Eigen::VectorXd & getDFS()        const { return dfs_; }
  const Eigen::VectorXd & getUsageCount() const { return usage_; }

 private:
  /// Ensure scratch buffers have shapes matching current (nens, nloc).
  void ensureBuffers(index_t nens, index_t nloc) const {
    if (Ya_buf_.rows() != nens || Ya_buf_.cols() != nloc) Ya_buf_.resize(nens, nloc);
    if (diagS_buf_.size() != nloc)                        diagS_buf_.resize(nloc);
  }

  // Cached global offsets for each block (recomputed when total nobs changes).
  const std::vector<std::size_t> & blockOffsetsCached(const Departures_ & omb) const {
    const std::size_t nblocks = omb.size();
    const std::size_t total   = omb.nobs();
    if (block_off_cache_.size() != nblocks || block_off_total_nobs_ != total) {
      block_off_cache_.resize(nblocks);
      std::size_t acc = 0;
      for (std::size_t b = 0; b < nblocks; ++b) {
        block_off_cache_[b] = acc;
        acc += static_cast<std::size_t>(omb[b].nobs());
      }
      block_off_total_nobs_ = total;
    }
    return block_off_cache_;
  }

  // Obtain obs types from YAML
  static std::vector<std::string> extractObsTypes(const eckit::Configuration & conf) {
    std::vector<eckit::LocalConfiguration> items = conf.getSubConfigurations();
    std::vector<std::string> labels;
    labels.reserve(items.size());
    for (const auto & sub : items) {
      if (!sub.has("obs space")) continue;
      eckit::LocalConfiguration obsSpace(sub, "obs space");
      if (obsSpace.has("obs type")) {
        labels.push_back(obsSpace.getString("obs type"));
      } else if (obsSpace.has("name")) {
        labels.push_back(obsSpace.getString("name"));
      } else {
        oops::Log::info() << "extractObsTypes: no 'obs type' or 'name' found in "
                          << "'obs space' in YAML; using 'Unknown'" << std::endl;
        labels.push_back("Unknown");
      }
    }
    return labels;
  }

 private:
  // Accumulators (global serial layout)
  Eigen::VectorXd dfs_;    // finalized global DFS after finalize()
  Eigen::VectorXd usage_;  // usage counts (for finalize())
  // Reusable temporaries
  mutable Eigen::MatrixXf Ya_buf_;
  mutable Eigen::VectorXd diagS_buf_;
  // Caches
  mutable std::vector<std::size_t> block_off_cache_;
  mutable std::size_t block_off_total_nobs_{std::numeric_limits<std::size_t>::max()};
};

}  // namespace oops
#endif  // OOPS_ASSIMILATION_DFSCALCULATOR_H_
