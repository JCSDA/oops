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
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace util {

// -----------------------------------------------------------------------------

void zeroFieldSet(atlas::FieldSet & fset) {
  oops::Log::trace() << "zeroFieldSet starting" << std::endl;

  for (auto & field : fset) {
    // Set data to zero
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
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

  oops::Log::trace() << "addFieldSets done" << std::endl;
}

// -----------------------------------------------------------------------------

void subtractFieldSets(atlas::FieldSet & fset,
                       const atlas::FieldSet & subFset) {
  oops::Log::trace() << "subtractFieldSets starting" << std::endl;

  // Loop over subtracted fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (auto & subField : subFset) {
    // Get field with the same name
    atlas::Field field = fset.field(subField.name());

    // Get data and sub
    if (field.rank() == 2 && subField.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto subView = atlas::array::make_view<double, 2>(subField);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) -= subView(jnode, jlevel);
        }
      }
    } else {
      ABORT("subFieldSets: wrong rank");
    }
  }

  oops::Log::trace() << "subFieldSets done" << std::endl;
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

  oops::Log::trace() << "multiplyFieldSets done" << std::endl;
}

// -----------------------------------------------------------------------------

double dotProductFieldsLocal(const atlas::Field & field1,
                             const atlas::Field & field2,
                             const bool & includeHalo) {
  oops::Log::trace() << "dotProductFieldsLocal starting" << std::endl;
  ASSERT(field1.name() == field2.name());
  // Compute dot product
  double dp = 0.0;
  if (field1.rank() == 2 && field2.rank() == 2) {
    // Check fields consistency
    ASSERT(field1.shape(0) == field2.shape(0));
    ASSERT(field1.shape(1) == field2.shape(1));

    // Add contributions
    auto view1 = atlas::array::make_view<double, 2>(field1);
    auto view2 = atlas::array::make_view<double, 2>(field2);
    if (field1.functionspace().type() == "Spectral") {
      const atlas::functionspace::Spectral fs(field1.functionspace());
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
              if (view1(jnode, jlevel) != util::missingValue<double>()
                && view2(jnode, jlevel) != util::missingValue<double>()) {
                dp += view1(jnode, jlevel)*view2(jnode, jlevel);
              }
            }
            ++jnode;

            // No imaginary part
            ++jnode;
          } else {
            // Real part
            for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
              if (view1(jnode, jlevel) != util::missingValue<double>()
                && view2(jnode, jlevel) != util::missingValue<double>()) {
                dp += 2.0*view1(jnode, jlevel)*view2(jnode, jlevel);
              }
            }
            ++jnode;

            // Imaginary part
            for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
              if (view1(jnode, jlevel) != util::missingValue<double>()
                && view2(jnode, jlevel) != util::missingValue<double>()) {
                dp += 2.0*view1(jnode, jlevel)*view2(jnode, jlevel);
              }
            }
            ++jnode;
          }
        }
      }
    } else {
      const auto ghostView = atlas::array::make_view<int, 1>(field1.functionspace().ghost());
      for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
        if (includeHalo || (ghostView(jnode) == 0)) {
          for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
            if (view1(jnode, jlevel) != util::missingValue<double>()
              && view2(jnode, jlevel) != util::missingValue<double>()) {
              dp += view1(jnode, jlevel)*view2(jnode, jlevel);
            }
          }
        }
      }
    }
  } else {
    ABORT("dotProductFieldsLocal: wrong rank");
  }

  oops::Log::trace() << "dotProductFieldsLocal done" << std::endl;
  // Return dot product
  return dp;
}

// -----------------------------------------------------------------------------

double dotProductFields(const atlas::Field & field1,
                        const atlas::Field & field2,
                        const eckit::mpi::Comm & comm,
                        const bool & includeHalo) {
  double dp = dotProductFieldsLocal(field1, field2, includeHalo);
  // Allreduce
  comm.allReduceInPlace(dp, eckit::mpi::sum());
  return dp;
}

// -----------------------------------------------------------------------------

double dotProductFieldSets(const atlas::FieldSet & fset1,
                           const atlas::FieldSet & fset2,
                           const std::vector<std::string> & vars,
                           const eckit::mpi::Comm & comm,
                           const bool & includeHalo) {
  oops::Log::trace() << "dotProductFieldSets starting" << std::endl;

  // Compute dot product
  double dp = 0.0;
  for (const auto & var : vars) {
    // Check fields presence
    if (fset1.has(var) && fset2.has(var)) {
      dp += dotProductFieldsLocal(fset1.field(var), fset2.field(var), includeHalo);
    }
  }
  // Allreduce
  comm.allReduceInPlace(dp, eckit::mpi::sum());

  oops::Log::trace() << "dotProductFieldSets done" << std::endl;
  // Return dot product
  return dp;
}

// -----------------------------------------------------------------------------

double normField(const atlas::Field & field,
                 const eckit::mpi::Comm & comm) {
  return std::sqrt(dotProductFields(field, field, comm, false));
}

// -----------------------------------------------------------------------------

double normFieldSet(const atlas::FieldSet & fset,
                    const std::vector<std::string> & vars,
                    const eckit::mpi::Comm & comm) {
  return std::sqrt(dotProductFieldSets(fset, fset, vars, comm, false));
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
          if (std::abs(divView(jnode, jlevel)) > 0.0) {
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

  oops::Log::trace() << "divideFieldSets done" << std::endl;
}

// -----------------------------------------------------------------------------

void divideFieldSets(atlas::FieldSet & fset,
                     const  atlas::FieldSet & divFset,
                     const  atlas::FieldSet & maskFset) {
  oops::Log::trace() << "divideFieldSets with mask starting" << std::endl;

  // Loop over divider fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (const auto & divField : divFset) {
    // Get field with the same name
    atlas::Field field = fset.field(divField.name());

    // Get mask field with the same name
    atlas::Field mask = maskFset.field(divField.name());

    // Get data and divide
    if (field.rank() == 2 && divField.rank() == 2 && mask.rank() == 2) {
      const auto divView = atlas::array::make_view<double, 2>(divField);
      const auto maskView = atlas::array::make_view<double, 2>(mask);
      auto view = atlas::array::make_view<double, 2>(field);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (std::abs(maskView(jnode, jlevel)) > 0.0) {
            if (std::abs(divView(jnode, jlevel)) > 0.0) {
              view(jnode, jlevel) /= divView(jnode, jlevel);
            } else {
                ABORT("divideFieldSets with mask: divide by zero");
            }
          }
        }
      }
    } else {
      ABORT("divideFieldSets with mask: wrong rank");
    }
  }

  oops::Log::trace() << "divideFieldSets with mask done" << std::endl;
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

  oops::Log::trace() << "sqrtFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------
// if a fieldset does not have the field with name defined by string variable
// "fldname", it will create one using the functionspace from field defined by
// variable string "template_name" and assign it a value of 0.0
void addZeroFieldToFieldSet(const std::string & fldname,
                            const std::string & template_name,
                            atlas::FieldSet & fset) {
  oops::Log::trace() << "addZeroFieldToFieldSet starting" << std::endl;

  if ((!fset.has(fldname)) && fset.has(template_name) &&
      (fset[template_name].rank() == 2)) {
    atlas::Field t = fset[template_name].functionspace().createField<double>(
      atlas::option::name(fldname) |
      atlas::option::levels(fset[template_name].levels()));
    t.haloExchange();
    atlas::array::make_view<double, 2>(t).assign(0.0);
    fset.add(t);
  }

  oops::Log::trace() << "addZeroFieldToFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util
