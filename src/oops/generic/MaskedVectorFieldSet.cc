// (C) Crown Copyright 2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#include "oops/generic/MaskedVectorFieldSet.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <string>

#include "atlas/array/MakeView.h"
#include "eckit/exception/Exceptions.h"

namespace oops {

namespace {

// Computes and stores where each (field, masked location, level) falls
// within a flat vector laid out as [field][masked location][level].
// Shared by fieldSetToVector and vectorToFieldSet below.
class MaskedVectorLayout {
 public:
  MaskedVectorLayout(const atlas::FieldSet& fieldSet,
                     const std::vector<bool>& locationMask)
      : locationIndices_{[&] {
          auto indices = std::vector<std::size_t>{};
          indices.reserve(
              std::count(locationMask.begin(), locationMask.end(), true));
          for (std::size_t locationIdx = 0; locationIdx < locationMask.size();
               ++locationIdx) {
            if (locationMask[locationIdx]) {
              indices.push_back(locationIdx);
            }
          }
          return indices;
        }()},
        variableDisplacements_{[&] {
          auto displacements = std::vector<std::size_t>{};
          displacements.reserve(fieldSet.size() + 1);
          displacements.push_back(0);
          std::transform_inclusive_scan(
              fieldSet.begin(), fieldSet.end(),
              std::back_inserter(displacements), std::plus<>{},
              [&](const atlas::Field& field) {
                return locationMask.size() *
                       static_cast<std::size_t>(field.shape(1));
              });
          return displacements;
        }()} {
    for (const auto& field : fieldSet) {
      if (static_cast<std::size_t>(field.shape(0)) != locationIndices_.size()) {
        throw eckit::BadValue(
            "Field location dimension does not match number of "
            "unmasked locations");
      }
    }
  }

  std::size_t totalSize() const { return variableDisplacements_.back(); }

  // Copy data element-by-element between a matching field in fieldSet and
  // dataVector for every field in this layout, using function to transfer a
  // single (field value, vector value) pair.
  template <typename FieldSetT, typename VectorT, typename Functor>
  void copyData(FieldSetT& fieldSet, VectorT& dataVector,
                const Functor& copier) const {
    for (atlas::idx_t fieldIdx = 0; fieldIdx < fieldSet.size(); ++fieldIdx) {
      auto fieldView = atlas::array::make_view<double, 2>(fieldSet[fieldIdx]);

      const auto numLocs = fieldView.shape(0);
      const auto numLevels = fieldView.shape(1);
      const auto dataBeginIdx = variableDisplacements_[fieldIdx];

      for (atlas::idx_t loc = 0; loc < numLocs; ++loc) {
        const auto locationBeginIdx =
            dataBeginIdx + locationIndices_[loc] * numLevels;
        for (atlas::idx_t lev = 0; lev < numLevels; ++lev) {
          copier(fieldView(loc, lev), dataVector[locationBeginIdx + lev]);
        }
      }
    }
  }

 private:
  std::vector<std::size_t> locationIndices_{};
  std::vector<std::size_t> variableDisplacements_{};
};

}  // namespace

void fieldSetToVector(const std::vector<bool>& locationMask,
                      const atlas::FieldSet& fieldSet,
                      std::vector<double>& dataVector) {
  const auto layout = MaskedVectorLayout(fieldSet, locationMask);
  dataVector.resize(layout.totalSize());

  layout.copyData(fieldSet, dataVector,
                  [&](double fieldValue, double& vectorValue) {
                    vectorValue = fieldValue;
                  });
}

void vectorToFieldSet(const std::vector<bool>& locationMask,
                      const std::vector<double>& dataVector,
                      atlas::FieldSet& fieldSet) {
  const auto layout = MaskedVectorLayout(fieldSet, locationMask);
  if (layout.totalSize() != dataVector.size()) {
    throw eckit::BadValue(
        "Data vector size does not match variable displacements.", Here());
  }

  layout.copyData(fieldSet, dataVector,
                  [&](double& fieldValue, double vectorValue) {
                    fieldValue = vectorValue;
                  });
}

}  // namespace oops
