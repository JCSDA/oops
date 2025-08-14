/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FieldSetOperations.h"

#include <cmath>

#include "atlas/array.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/for_each.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/reduction.h"

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
      throw eckit::Exception("zeroFieldSet: wrong rank", Here());
    }
  }

  fset.set_dirty(false);

  oops::Log::trace() << "zeroFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void addFieldSets(atlas::FieldSet & fset,
                  const atlas::FieldSet & addFset) {
  oops::Log::trace() << "addFieldSets starting" << std::endl;

  // Loop over additive fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (auto & addField : addFset) {
    atlas::Field field = fset.field(addField.name());
    util::for_each_value(
      util::IndexRange::include_halo,  // atlas 0.43 will enable excluding for all FunctionSpaces
      [](const double rhs, double & lhs) { lhs += rhs; },
      addField,
      field);

    // If either term in the sum is out-of-date, then the result will be out-of-date
    field.set_dirty(field.dirty() || addField.dirty());
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
    atlas::Field field = fset.field(subField.name());
    util::for_each_value(
      util::IndexRange::include_halo,  // atlas 0.43 will enable excluding for all FunctionSpaces
      [](const double rhs, double & lhs) { lhs -= rhs; },
      subField,
      field);

    // If either term in the subtraction is out-of-date, then the result will be out-of-date
    field.set_dirty(field.dirty() || subField.dirty());
  }

  oops::Log::trace() << "subFieldSets done" << std::endl;
}

// -----------------------------------------------------------------------------

void multiplyFieldSet(atlas::FieldSet & fset,
                      const double mul) {
  oops::Log::trace() << "multiplyFieldSet starting" << std::endl;

  // Loop over fields
  for (auto & field : fset) {
    util::for_each_value(
      util::IndexRange::include_halo,  // atlas 0.43 will enable excluding for all FunctionSpaces
      [=](double & val) { val *= mul; },
      field);
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
    atlas::Field field = fset.field(mulField.name());
    util::for_each_value(
      util::IndexRange::include_halo,  // atlas 0.43 will enable excluding for all FunctionSpaces
      [](const double rhs, double & lhs) { lhs *= rhs; },
      mulField,
      field);

    // If either term in the product is out-of-date, then the result will be out-of-date
    field.set_dirty(field.dirty() || mulField.dirty());
  }

  oops::Log::trace() << "multiplyFieldSets done" << std::endl;
}

// -----------------------------------------------------------------------------

double dotProductFields(const atlas::Field & field1,
                        const atlas::Field & field2,
                        const eckit::mpi::Comm & comm) {
  // Compute task-local dot product, then sum over tasks
  double dp = util::dot_product_on_task(field1, field2);

  comm.allReduceInPlace(dp, eckit::mpi::sum());

  return dp;
}

// -----------------------------------------------------------------------------

double dotProductFieldSets(const atlas::FieldSet & fset1,
                           const atlas::FieldSet & fset2,
                           const std::vector<std::string> & vars,
                           const eckit::mpi::Comm & comm) {
  oops::Log::trace() << "dotProductFieldSets starting" << std::endl;

  // Compute task-local dot product for Fields in both FieldSets, then sum over tasks
  double dp = 0.0;
  for (const auto & var : vars) {
    if (fset1.has(var) && fset2.has(var)) {
      dp += util::dot_product_on_task(fset1.field(var), fset2.field(var));
    }
  }

  comm.allReduceInPlace(dp, eckit::mpi::sum());

  oops::Log::trace() << "dotProductFieldSets done" << std::endl;
  return dp;
}

// -----------------------------------------------------------------------------

double normField(const atlas::Field & field,
                 const eckit::mpi::Comm & comm) {
  return std::sqrt(dotProductFields(field, field, comm));
}

// -----------------------------------------------------------------------------

double normFieldSet(const atlas::FieldSet & fset,
                    const std::vector<std::string> & vars,
                    const eckit::mpi::Comm & comm) {
  return std::sqrt(dotProductFieldSets(fset, fset, vars, comm));
}

// -----------------------------------------------------------------------------

void divideFieldSets(atlas::FieldSet & fset,
                     const  atlas::FieldSet & divFset) {
  oops::Log::trace() << "divideFieldSets starting" << std::endl;

  // Loop over divider fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (const auto & divField : divFset) {
    bool found_div_by_zero = false;
    atlas::Field field = fset.field(divField.name());
    util::for_each_value(
      util::IndexRange::include_halo,  // atlas 0.43 will enable excluding for all FunctionSpaces
      [&](const double rhs, double & lhs) {
        if (std::abs(rhs) > 0.0) {
          lhs /= rhs;
        } else if (std::abs(lhs) > 0.0) {
          // If the numerator is 0, then it's ok for the denominator to be 0; this is probably
          // a case of 0/0 in the halo, and we opt to return 0 (i.e., no change to the field).
          // If the numerator is finite (= this else branch), this is a divide-by-zero error:
          found_div_by_zero = true;
        }
      },
      divField,
      field);

    if (found_div_by_zero) {
      throw eckit::Exception("divideFieldSets: divide by zero for field " + divField.name(),
                             Here());
    }

    // If either term in the division is out-of-date, then the result will be out-of-date
    field.set_dirty(field.dirty() || divField.dirty());
  }

  oops::Log::trace() << "divideFieldSets done" << std::endl;
}

// -----------------------------------------------------------------------------

void divideFieldSets(atlas::FieldSet & fset,
                     const atlas::FieldSet & divFset,
                     const atlas::FieldSet & maskFset) {
  oops::Log::trace() << "divideFieldSets with mask starting" << std::endl;

  // Loop over divider fields. The RHS FieldSet may contain only a subset of Fields from the
  // input/output FieldSet. If this is the case, no work is done for fields present only in the LHS.
  for (const auto & divField : divFset) {
    bool found_div_by_zero = false;
    atlas::Field field = fset.field(divField.name());
    atlas::Field mask = maskFset.field(divField.name());
    util::for_each_value(
      util::IndexRange::include_halo,  // atlas 0.43 will enable excluding for all FunctionSpaces
      [&](const double rhs, const double mask, double & lhs) {
        if (std::abs(mask) > 0.0) {
          if (std::abs(rhs) > 0.0) {
            lhs /= rhs;
          } else if (std::abs(lhs) > 0.0) {
            found_div_by_zero = true;
          }
        }
      },
      divField,
      mask,
      field);

    if (found_div_by_zero) {
      throw eckit::Exception("divideFieldSets: divide by zero for field " + divField.name(),
                             Here());
    }

    // If either term in the division is out-of-date, then the result will be out-of-date
    field.set_dirty(field.dirty() || divField.dirty());
  }

  oops::Log::trace() << "divideFieldSets with mask done" << std::endl;
}

// -----------------------------------------------------------------------------

void sqrtFieldSet(atlas::FieldSet & fset) {
  oops::Log::trace() << "sqrtFieldSet starting" << std::endl;

  // Loop over fields
  for (auto field : fset) {
    bool found_negative_sqrt = false;
    util::for_each_value(
      [&](double & val) {
        if (val >= 0.0) {
          val = std::sqrt(val);
        } else {
          found_negative_sqrt = true;
        }
      },
      field);

    if (found_negative_sqrt) {
      throw eckit::Exception("sqrtFieldSet: negative square root for field " + field.name(),
                             Here());
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
      atlas::option::levels(fset[template_name].shape(1)));
    t.haloExchange();
    atlas::array::make_view<double, 2>(t).assign(0.0);
    fset.add(t);
  }

  oops::Log::trace() << "addZeroFieldToFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util
