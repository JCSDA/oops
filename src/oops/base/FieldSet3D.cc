/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/FieldSet3D.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

FieldSet3D initFieldSet3D(const FieldSet3D & other) {
  FieldSet3D copy(other.validTime(), other.commGeom());
  for (const auto & otherField : other.fieldSet()) {
    // Create Field
    atlas::Field field = otherField.functionspace().createField<double>(
      atlas::option::name(otherField.name()) | atlas::option::levels(otherField.levels()));
    // Copy metadata
    field.metadata() = otherField.metadata();
    // Add field
    copy.fieldSet().add(field);
  }
  return copy;
}

// -----------------------------------------------------------------------------

FieldSet3D::FieldSet3D(const atlas::FieldSet & fset,
                       const util::DateTime & validTime,
                       const eckit::mpi::Comm & comm):
  fset_(fset), validTime_(validTime), vars_(currentVariables()),
  comm_(comm)
{}

// -----------------------------------------------------------------------------

FieldSet3D::FieldSet3D(const util::DateTime & validTime,
                       const eckit::mpi::Comm & comm):
  fset_(), validTime_(validTime), vars_(), comm_(comm)
{}
// -----------------------------------------------------------------------------

oops::Variables FieldSet3D::currentVariables() const {
  oops::Variables vars;
  for (const auto & field : fset_) {
    vars.push_back(field.name());
    vars.addMetaData(field.name(), "levels", field.levels());
  }
  return vars;
}

// -----------------------------------------------------------------------------

const oops::Variables & FieldSet3D::variables() const {
  vars_ = currentVariables();
  return vars_;
}

// -----------------------------------------------------------------------------

void FieldSet3D::zero() {
  util::zeroFieldSet(fset_);
}

// -----------------------------------------------------------------------------

FieldSet3D & FieldSet3D::operator+=(const FieldSet3D & other) {
  util::addFieldSets(fset_, other.fset_);
  return *this;
}

// -----------------------------------------------------------------------------

double FieldSet3D::dot_product_with(const FieldSet3D & other, const oops::Variables & vars) const {
  return util::dotProductFieldSets(fset_, other.fieldSet(), vars.variables(), comm_);
}

// -----------------------------------------------------------------------------

void FieldSet3D::print(std::ostream & os) const {
  os << std::endl << "Valid time: " << validTime() << std::endl;
  for (const auto & field : fset_) {
    os << " Field " << field.name()
       << " norm: " << util::normField(field, comm_) << std::endl;
  }
}

// -----------------------------------------------------------------------------

size_t FieldSet3D::serialSize() const {
  // size of valid time + number of variables
  size_t fset_size = validTime_.serialSize() + 1;
  for (const auto & field : fset_) {
    assert(field.rank() == 2);
    // size of field + dimension sizes (2) + variable name hash (1)
    fset_size += field.shape(0) * field.shape(1) + 3;
  }
  return fset_size;
}

// -----------------------------------------------------------------------------

void FieldSet3D::serialize(std::vector<double> & vect)  const {
  // Allocate space
  const size_t fset_size = this->serialSize();
  vect.reserve(vect.size() + fset_size);

  // serialize valid time and number of variables
  validTime_.serialize(vect);
  vect.push_back(fset_.size());

  // Serialize the fields, including variable name hashes and sizes
  for (const auto & field : fset_) {
    const size_t varname_hash = std::hash<std::string>{}(field.name());
    vect.push_back(reinterpret_cast<const double &>(varname_hash));
    vect.push_back(field.shape(0));
    vect.push_back(field.shape(1));
    size_t index = vect.size();
    vect.resize(vect.size() + field.shape(0) * field.shape(1));
    const auto view = atlas::array::make_view<double, 2>(field);
    for (int jnode = 0; jnode < field.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
        vect[index++] = view(jnode, jlevel);
      }
    }
  }
}

// -----------------------------------------------------------------------------

void FieldSet3D::deserialize(const std::vector<double> & vect, size_t & index) {
  util::DateTime other_time;
  other_time.deserialize(vect, index);
  if (other_time != validTime_) {
    // All current use cases for this method are needed for fieldsets at different
    // times to handle 4D aspects in covariances: issue a warning that the dates are
    // different (may be useful for future use cases) but proceed.
    oops::Log::warning() << "FieldSet3D::deserialize: valid times are different." << std::endl;
    oops::Log::warning() << "This FieldSet3D valid time: " << validTime_ << std::endl;
    oops::Log::warning() << "Valid time of the serialized FieldSet3D: " << other_time << std::endl;
  }
  const int other_nvars = vect[index++];
  if (other_nvars != fset_.size()) {
    oops::Log::error() << "This FieldSet3D number of variables: " << fset_.size() << std::endl;
    oops::Log::error() << "Serialized FieldSet3D number of variables: " << other_nvars << std::endl;
    throw eckit::BadParameter("FieldSet3D::deserialize failed: different number of variables",
                              Here());
  }
  // Deserialize the fields
  for (auto & field : fset_) {
    const size_t other_varname_hash = reinterpret_cast<const size_t &>(vect[index++]);
    if (other_varname_hash != std::hash<std::string>{}(field.name())) {
      oops::Log::error() << "This FieldSet3D variable " << field.name() << " does not match "
                         << "corresponding serialized FieldSet3D variable." << std::endl;
      throw eckit::BadParameter("FieldSet3D::deserialize failed: different variable name", Here());
    }
    const int other_field_shape_0 = vect[index++];
    const int other_field_shape_1 = vect[index++];
    if (other_field_shape_0 != field.shape(0) || other_field_shape_1 != field.shape(1)) {
      oops::Log::error() << "This FieldSet3D variable " << field.name() << " sizes: "
                         << field.shape(0) << ", " << field.shape(1) << std::endl;
      oops::Log::error() << "Serialized FieldSet3D variable " << field.name() << " sizes: "
                         << other_field_shape_0 << ", " << other_field_shape_1 << std::endl;
      throw eckit::BadParameter("FieldSet3D::deserialize failed: different field shapes", Here());
    }
    auto view = atlas::array::make_view<double, 2>(field);
    for (int jnode = 0; jnode < field.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
        view(jnode, jlevel) = vect[index++];
      }
    }
  }
}

}  // namespace oops
