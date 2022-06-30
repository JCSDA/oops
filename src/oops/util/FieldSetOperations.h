/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_FIELDSETOPERATIONS_H_
#define OOPS_UTIL_FIELDSETOPERATIONS_H_

#include "atlas/field.h"

#include "oops/util/Logger.h"

namespace util {

// -----------------------------------------------------------------------------

void FieldSetMultiply(atlas::FieldSet & fset, const atlas::FieldSet & mulFset) {
  oops::Log::trace() << "FieldSetMultiply starting" << std::endl;
  // Loop over multiplier fields
  for (const auto mulField : mulFset) {
    // Get increment field with the same name
    atlas::Field field = fset.field(mulField.name());

    // Get data and multiply
    if (mulField.rank() == 1) {
      const auto mulView = atlas::array::make_view<double, 1>(mulField);
      auto view = atlas::array::make_view<double, 1>(field);
      for (int j0 = 0; j0 < mulField.shape(0); ++j0) {
        view(j0) *= mulView(j0);
      }
    }
    if (mulField.rank() == 2) {
      const auto mulView = atlas::array::make_view<double, 2>(mulField);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int j0 = 0; j0 < mulField.shape(0); ++j0) {
        for (int j1 = 0; j1 < mulField.shape(1); ++j1) {
          view(j0, j1) *= mulView(j0, j1);
        }
      }
    }
    if (mulField.rank() == 3) {
      const auto mulView = atlas::array::make_view<double, 3>(mulField);
      auto view = atlas::array::make_view<double, 3>(field);
      for (int j0 = 0; j0 < mulField.shape(0); ++j0) {
        for (int j1 = 0; j1 < mulField.shape(1); ++j1) {
          for (int j2 = 0; j2 < mulField.shape(2); ++j2) {
            view(j0, j1, j2) *= mulView(j0, j1, j2);
          }
        }
      }
    }
  }
  oops::Log::trace() << "FieldSetMultiply starting" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldSetDivide(atlas::FieldSet & fset, const  atlas::FieldSet & divFset) {
  oops::Log::trace() << "FieldSetDivide starting" << std::endl;
  // Loop over divider fields
  for (const auto divField : divFset) {
    // Get increment field with the same name
    atlas::Field field = fset.field(divField.name());

    // Get data and divide
    if (divField.rank() == 1) {
      const auto divView = atlas::array::make_view<double, 1>(divField);
      auto view = atlas::array::make_view<double, 1>(field);
      for (int j0 = 0; j0 < divField.shape(0); ++j0) {
        view(j0) /= divView(j0);
      }
    }
    if (divField.rank() == 2) {
      const auto divView = atlas::array::make_view<double, 2>(divField);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int j0 = 0; j0 < divField.shape(0); ++j0) {
        for (int j1 = 0; j1 < divField.shape(1); ++j1) {
          view(j0, j1) /= divView(j0, j1);
        }
      }
    }
    if (divField.rank() == 3) {
      const auto divView = atlas::array::make_view<double, 3>(divField);
      auto view = atlas::array::make_view<double, 3>(field);
      for (int j0 = 0; j0 < divField.shape(0); ++j0) {
        for (int j1 = 0; j1 < divField.shape(1); ++j1) {
          for (int j2 = 0; j2 < divField.shape(2); ++j2) {
            view(j0, j1, j2) /= divView(j0, j1, j2);
          }
        }
      }
    }
  }
  oops::Log::trace() << "FieldSetDivide starting" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldSetSqrt(atlas::FieldSet & fset) {
  oops::Log::trace() << "FieldSetSqrt starting" << std::endl;
  // Loop over fields
  for (auto field : fset) {
    // Get data and take square-root
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      for (int j0 = 0; j0 < field.shape(0); ++j0) {
        ASSERT(view(j0) >= 0.0);
        view(j0) = std::sqrt(view(j0));
      }
    }
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (int j0 = 0; j0 < field.shape(0); ++j0) {
        for (int j1 = 0; j1 < field.shape(1); ++j1) {
          ASSERT(view(j0, j1) >= 0.0);
          view(j0, j1) = std::sqrt(view(j0, j1));
        }
      }
    }
    if (field.rank() == 3) {
      auto view = atlas::array::make_view<double, 3>(field);
      for (int j0 = 0; j0 < field.shape(0); ++j0) {
        for (int j1 = 0; j1 < field.shape(1); ++j1) {
          for (int j2 = 0; j2 < field.shape(2); ++j2) {
            ASSERT(view(j0, j1, j2) >= 0.0);
            view(j0, j1, j2) = std::sqrt(view(j0, j1, j2));
          }
        }
      }
    }
  }
  oops::Log::trace() << "FieldSetSqrt starting" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_FIELDSETOPERATIONS_H_
