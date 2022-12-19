/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FieldSetOperations.h"

#include <cmath>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

namespace util {

// -----------------------------------------------------------------------------
atlas::FieldSet createRandomFieldSet(const oops::GeometryData & geometryData,
                                     const std::vector<size_t> & variableSizes,
                                     const oops::Variables & vars) {
  oops::Log::trace() << "createRandomFieldSet starting" << std::endl;

  // Get ghost points
  atlas::Field ghost;
  if (geometryData.functionSpace().type() == "Spectral") {
    ghost = geometryData.functionSpace().createField<int>(atlas::option::name("ghost"));
    auto ghostView = atlas::array::make_view<int, 1>(ghost);
    for (atlas::idx_t jnode = 0; jnode < ghost.shape(0); ++jnode) {
      ghostView(jnode) = 0;
    }
  } else {
    ghost = geometryData.functionSpace().ghost();
  }
  auto ghostView = atlas::array::make_view<int, 1>(ghost);

  // Create FieldSet
  atlas::FieldSet fset;

  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    // Create field
    atlas::Field field = geometryData.functionSpace().createField<double>(
      atlas::option::name(vars.variables()[jvar]) | atlas::option::levels(variableSizes[jvar]));

    // Get field owned size
    size_t n = 0;
    if (field.rank() == 2) {
      if (geometryData.functionSpace().type() == "Spectral") {
        atlas::functionspace::Spectral fs(geometryData.functionSpace());
        const atlas::idx_t N = fs.truncation();
        const auto zonal_wavenumbers = fs.zonal_wavenumbers();
        const atlas::idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();
        for (int jm=0; jm < nb_zonal_wavenumbers; ++jm) {
          const atlas::idx_t m1 = zonal_wavenumbers(jm);
          for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
            if (m1 == 0) {
              n += field.shape(1);
            } else {
              n += 2*field.shape(1);
            }
          }
        }
      } else {
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          if (ghostView(jnode) == 0) n += field.shape(1);
        }
      }
    }

    // Generate random vector
    util::NormalDistribution<double> rand_vec(n, 0.0, 1.0, 1);

    // Populate with random numbers
    n = 0;
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      if (geometryData.functionSpace().type() == "Spectral") {
        atlas::functionspace::Spectral fs(geometryData.functionSpace());
        const atlas::idx_t N = fs.truncation();
        const auto zonal_wavenumbers = fs.zonal_wavenumbers();
        const atlas::idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();
        int jnode = 0;
        for (int jm=0; jm < nb_zonal_wavenumbers; ++jm) {
          const atlas::idx_t m1 = zonal_wavenumbers(jm);
          for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
            if (m1 == 0) {
              // Real part only
              for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                view(jnode, jlevel) = rand_vec[n];
                ++n;
              }
              ++jnode;

              // No imaginary part
              for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                view(jnode, jlevel) = 0.0;
              }
              ++jnode;
            } else {
              // Real part
              for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                view(jnode, jlevel) = rand_vec[n];
                ++n;
              }
              ++jnode;

              // Imaginary part
              for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                view(jnode, jlevel) = rand_vec[n];
                ++n;
              }
              ++jnode;
            }
          }
        }
      } else {
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            if (ghostView(jnode) == 0) {
              view(jnode, jlevel) = rand_vec[n];
              ++n;
            } else {
              view(jnode, jlevel) = 0.0;
            }
          }
        }
      }
    }

    // Add field
    fset.add(field);
  }

  if (geometryData.functionSpace().type() != "Spectral") {
    // Halo exchange
    fset.haloExchange();
  }

  // Return FieldSet
  return fset;
}

// -----------------------------------------------------------------------------

atlas::FieldSet copyFieldSet(const atlas::FieldSet & otherFset) {
  oops::Log::trace() << "copyFieldSet starting" << std::endl;

  // Create FieldSet
  atlas::FieldSet fset;

  for (const auto & otherField : otherFset) {
    // Create Field
    atlas::Field field = otherField.functionspace().createField<double>(
      atlas::option::name(otherField.name()) | atlas::option::levels(otherField.levels()));

    // Copy data
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto otherView = atlas::array::make_view<double, 2>(otherField);
      for (atlas::idx_t jnode = 0; jnode < otherField.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < otherField.shape(1); ++jlevel) {
            view(jnode, jlevel) = otherView(jnode, jlevel);
        }
      }
    } else {
      ABORT("copyFieldSet: wrong rank");
    }

    // Add field
    fset.add(field);
  }

  // Return FieldSet
  return fset;
}

// -----------------------------------------------------------------------------

void removeFieldsFromFieldSet(atlas::FieldSet & fset,
                              const oops::Variables & vars) {
  oops::Log::trace() << "removeFieldsFromFieldSet starting" << std::endl;

  // Create FieldSet
  atlas::FieldSet fsetTmp;

  for (const auto & field : fset) {
    if (!vars.has(field.name())) {
      // Add field
      fsetTmp.add(field);
    }
  }

  // Replace FieldSet
  fset = fsetTmp;

  oops::Log::trace() << "removeFieldsFromFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void zeroFieldSet(atlas::FieldSet & fset) {
  oops::Log::trace() << "zeroFieldSet starting" << std::endl;

  for (auto & field : fset) {
    // Set data to zero
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) = 0.0;
        }
      }
    } else {
      ABORT("zeroFieldSet: wrong rank");
    }
  }

  oops::Log::trace() << "zeroFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void addFieldSets(atlas::FieldSet & fset,
                  const atlas::FieldSet & addFset) {
  oops::Log::trace() << "addFieldSets starting" << std::endl;

  // Loop over additive fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (auto & addField : addFset) {
    // Get field with the same name
    atlas::Field field = fset.field(addField.name());

    // Get data and add
    if (field.rank() == 2 && addField.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto addView = atlas::array::make_view<double, 2>(addField);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) += addView(jnode, jlevel);
        }
      }
    } else {
      ABORT("addFieldSets: wrong rank");
    }
  }

  oops::Log::trace() << "addFieldSets starting" << std::endl;
}

// -----------------------------------------------------------------------------

void multiplyFieldSet(atlas::FieldSet & fset,
                      const double & mul) {
  oops::Log::trace() << "multiplyFieldSet starting" << std::endl;

  // Loop over fields
  for (auto & field : fset) {
    // Get data and multiply
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= mul;
        }
      }
    } else {
      ABORT("multiplyFieldSet: wrong rank");
    }
  }
  oops::Log::trace() << "multiplyFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void multiplyFieldSets(atlas::FieldSet & fset,
                       const atlas::FieldSet & mulFset) {
  oops::Log::trace() << "multiplyFieldSets starting" << std::endl;

  // Loop over multiplier fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (const auto & mulField : mulFset) {
    // Get field with the same name
    atlas::Field field = fset.field(mulField.name());

    // Get data and multiply
    if (field.rank() == 2 && mulField.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto mulView = atlas::array::make_view<double, 2>(mulField);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= mulView(jnode, jlevel);
        }
      }
    } else {
      ABORT("multiplyFieldSets: wrong rank");
    }
  }

  oops::Log::trace() << "multiplyFieldSets starting" << std::endl;
}

// -----------------------------------------------------------------------------

double dotProductFieldSets(const atlas::FieldSet & fset1,
                           const atlas::FieldSet & fset2,
                           const oops::Variables & activeVars,
                           const eckit::mpi::Comm & comm) {
  oops::Log::trace() << "dotProductFieldSets starting" << std::endl;

  // Compute dot product
  double dp = 0.0;
  for (const auto & var : activeVars.variables()) {
    // Check fields presence
    if (fset1.has(var) && fset2.has(var)) {
      // Get fields
      const auto field1 = fset1.field(var);
      const auto field2 = fset2.field(var);

      if (field1.rank() == 2 && field2.rank() == 2) {
      // Check fields consistency
        ASSERT(field1.shape(0) == field2.shape(0));
        ASSERT(field1.shape(1) == field2.shape(1));

        // Add contributions
        auto view1 = atlas::array::make_view<double, 2>(field1);
        auto view2 = atlas::array::make_view<double, 2>(field2);
        if (field1.functionspace().type() == "Spectral") {
          atlas::functionspace::Spectral fs(field1.functionspace());
          const atlas::idx_t N = fs.truncation();
          const auto zonal_wavenumbers = fs.zonal_wavenumbers();
          const atlas::idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();
          int jnode = 0;
          for (int jm=0; jm < nb_zonal_wavenumbers; ++jm) {
            const atlas::idx_t m1 = zonal_wavenumbers(jm);
            for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
              if (m1 == 0) {
                // Real part only
                for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
                  dp += view1(jnode, jlevel)*view2(jnode, jlevel);
                }
                ++jnode;

                // No imaginary part
                ++jnode;
              } else {
                // Real part
                for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
                  dp += 2.0*view1(jnode, jlevel)*view2(jnode, jlevel);
                }
                ++jnode;

                // Imaginary part
                for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
                  dp += 2.0*view1(jnode, jlevel)*view2(jnode, jlevel);
                }
                ++jnode;
              }
            }
          }
        } else {
          for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
            for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
              dp += view1(jnode, jlevel)*view2(jnode, jlevel);
            }
          }
        }
      } else {
        ABORT("dotProductFieldSets: wrong rank");
      }
    }
  }

  // Allreduce
  comm.allReduceInPlace(dp, eckit::mpi::sum());

  // Return dot product
  return dp;

  oops::Log::trace() << "dotProductFieldSets done" << std::endl;
}

// -----------------------------------------------------------------------------

void divideFieldSets(atlas::FieldSet & fset,
                     const  atlas::FieldSet & divFset) {
  oops::Log::trace() << "divideFieldSets starting" << std::endl;

  // Loop over divider fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (const auto & divField : divFset) {
    // Get field with the same name
    atlas::Field field = fset.field(divField.name());

    // Get data and divide
    if (field.rank() == 2 && divField.rank() == 2) {
      const auto divView = atlas::array::make_view<double, 2>(divField);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (abs(divView(jnode, jlevel)) > 0.0) {
            view(jnode, jlevel) /= divView(jnode, jlevel);
          } else {
            ABORT("divideFieldSets: divide by zero");
          }
        }
      }
    } else {
      ABORT("divideFieldSets: wrong rank");
    }
  }

  oops::Log::trace() << "divideFieldSets starting" << std::endl;
}

// -----------------------------------------------------------------------------

void sqrtFieldSet(atlas::FieldSet & fset) {
  oops::Log::trace() << "sqrtFieldSet starting" << std::endl;

  // Loop over fields
  for (auto field : fset) {
    // Get data and take square-root
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          ASSERT(view(jnode, jlevel) >= 0.0);
          view(jnode, jlevel) = std::sqrt(view(jnode, jlevel));
        }
      }
    } else {
      ABORT("sqrtFieldSet: wrong rank");
    }
  }

  oops::Log::trace() << "sqrtFieldSet starting" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util

