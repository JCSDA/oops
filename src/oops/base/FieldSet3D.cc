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

#include "oops/base/Variable.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

FieldSet3D initFieldSet3D(const FieldSet3D & other) {
  FieldSet3D copy(other.validTime(), other.commGeom());
  copy.allocateOnly(other.fieldSet());
  return copy;
}

// -----------------------------------------------------------------------------

FieldSet3D initFieldSet3D(const util::DateTime & validTime,
                          const eckit::mpi::Comm & comm,
                          const atlas::FieldSet & other) {
  FieldSet3D copy(validTime, comm);
  copy.allocateOnly(other);
  return copy;
}

// -----------------------------------------------------------------------------

FieldSet3D copyFieldSet3D(const FieldSet3D & other,
                          const bool & shallow) {
  FieldSet3D copy(other.validTime(), other.commGeom());
  if (shallow) {
    copy.shallowCopy(other.fieldSet());
  } else {
    copy.deepCopy(other.fieldSet());
  }
  return copy;
}


// -----------------------------------------------------------------------------

FieldSet3D copyFieldSet3D(const util::DateTime & validTime,
                          const eckit::mpi::Comm & comm,
                          const atlas::FieldSet & other,
                          const bool & shallow) {
  FieldSet3D copy(validTime, comm);
  if (shallow) {
    copy.shallowCopy(other);
  } else {
    copy.deepCopy(other);
  }
  return copy;
}

// -----------------------------------------------------------------------------

FieldSet3D randomFieldSet3D(const util::DateTime & validTime,
                            const eckit::mpi::Comm & comm,
                            const atlas::FunctionSpace & fspace,
                            const Variables & vars) {
  FieldSet3D random(validTime, comm);
  random.randomInit(fspace, vars);
  return random;
}

// -----------------------------------------------------------------------------

FieldSet3D::FieldSet3D(const util::DateTime & validTime,
                       const eckit::mpi::Comm & comm):
  fset_(), validTime_(validTime), comm_(comm), name_()
{}

// -----------------------------------------------------------------------------

FieldSet3D::FieldSet3D(const FieldSet3D & other):
  fset_(), validTime_(other.validTime_), comm_(other.comm_), name_(other.name_) {
  fset_ = other.fset_.clone();
}

// -----------------------------------------------------------------------------

void FieldSet3D::allocateOnly(const atlas::FieldSet & otherFset) {
  fset_.clear();
  for (const auto & otherField : otherFset) {
    // Create field
    atlas::Field field = otherField.functionspace().createField<double>(
      atlas::option::name(otherField.name()) | atlas::option::levels(otherField.shape(1)));

    // Copy metadata
    field.metadata() = otherField.metadata();

    // Add field
    fset_.add(field);
  }
}

// -----------------------------------------------------------------------------

void FieldSet3D::allocateOnly(const FieldSet3D & otherFset) {
  ASSERT(validTime_ == otherFset.validTime_);
  this->allocateOnly(otherFset.fset_);
}

// -----------------------------------------------------------------------------

void FieldSet3D::deepCopy(const atlas::FieldSet & otherFset) {
  fset_ = otherFset.clone();
}

// -----------------------------------------------------------------------------

void FieldSet3D::deepCopy(const FieldSet3D & otherFset) {
  ASSERT(validTime_ == otherFset.validTime_);
  this->deepCopy(otherFset.fset_);
}

// -----------------------------------------------------------------------------

void FieldSet3D::shallowCopy(const atlas::FieldSet & otherFset) {
  fset_.clear();
  for (const auto & otherField : otherFset) {
    fset_.add(otherField);
  }
}

// -----------------------------------------------------------------------------

void FieldSet3D::shallowCopy(const FieldSet3D & otherFset) {
  ASSERT(validTime_ == otherFset.validTime_);
  this->shallowCopy(otherFset.fset_);
}

// -----------------------------------------------------------------------------

void FieldSet3D::init(const atlas::FunctionSpace & fspace,
                      const Variables & vars) {
  fset_ = util::createFieldSet(fspace, vars);
}

// -----------------------------------------------------------------------------

void FieldSet3D::init(const atlas::FunctionSpace & fspace,
                      const Variables & vars,
                      const double & value) {
  fset_ = util::createFieldSet(fspace, vars, value);
}

// -----------------------------------------------------------------------------

void FieldSet3D::randomInit(const atlas::FunctionSpace & fspace,
                            const Variables & vars) {
  fset_ = util::createRandomFieldSet(comm_, fspace, vars);
}

// -----------------------------------------------------------------------------

Variables FieldSet3D::currentVariables() const {
  Variables vars;
  for (const auto & field : fset_) {
    eckit::LocalConfiguration conf;
    conf.set("levels", field.shape(1));
    vars.push_back(Variable(field.name(), conf));
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

FieldSet3D & FieldSet3D::operator-=(const FieldSet3D & other) {
  util::subtractFieldSets(fset_, other.fset_);
  return *this;
}

// -----------------------------------------------------------------------------

FieldSet3D & FieldSet3D::operator*=(const FieldSet3D & other) {
  util::multiplyFieldSets(fset_, other.fset_);
  return *this;
}

// -----------------------------------------------------------------------------

FieldSet3D & FieldSet3D::operator*=(const double & zz) {
  util::multiplyFieldSet(fset_, zz);
  return *this;
}

// -----------------------------------------------------------------------------

FieldSet3D & FieldSet3D::operator/=(const FieldSet3D & other) {
  util::divideFieldSets(fset_, other.fset_);
  return *this;
}

// -----------------------------------------------------------------------------

double FieldSet3D::dot_product_with(const FieldSet3D & other, const Variables & vars) const {
  return util::dotProductFieldSets(fset_, other.fieldSet(), vars.variables(), comm_);
}

// -----------------------------------------------------------------------------

double FieldSet3D::norm(const Variables & vars) const {
  return util::normFieldSet(fset_, vars.variables(), comm_);
}

// -----------------------------------------------------------------------------

void FieldSet3D::sqrt() {
  util::sqrtFieldSet(fset_);
}

// -----------------------------------------------------------------------------

void FieldSet3D::read(const atlas::FunctionSpace & fspace,
                      const Variables & vars,
                      const eckit::LocalConfiguration & conf) {
  util::readFieldSet(comm_, fspace, vars, conf, fset_);
}

// -----------------------------------------------------------------------------

void FieldSet3D::write(const eckit::LocalConfiguration & conf) const {
  util::writeFieldSet(comm_, conf, fset_);
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

  static_assert(sizeof(double) == sizeof(size_t));
  // Serialize the fields, including variable name hashes and sizes
  for (const auto & field : fset_) {
    const size_t varname_hash = std::hash<std::string>{}(field.name());
    const double* varname_hash_asdouble = reinterpret_cast<const double*>(&varname_hash);
    vect.push_back(*varname_hash_asdouble);
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
    Log::warning() << "FieldSet3D::deserialize: valid times are different." << std::endl;
    Log::warning() << "This FieldSet3D valid time: " << validTime_ << std::endl;
    Log::warning() << "Valid time of the serialized FieldSet3D: " << other_time << std::endl;
  }
  const int other_nvars = vect[index++];
  if (other_nvars != fset_.size()) {
    Log::error() << "This FieldSet3D number of variables: " << fset_.size() << std::endl;
    Log::error() << "Serialized FieldSet3D number of variables: " << other_nvars << std::endl;
    throw eckit::BadParameter("FieldSet3D::deserialize failed: different number of variables",
                              Here());
  }
  // Deserialize the fields
  static_assert(sizeof(double) == sizeof(size_t));
  for (auto & field : fset_) {
    const size_t* other_varname_hash = reinterpret_cast<const size_t*>(&vect[index++]);
    if (*other_varname_hash != std::hash<std::string>{}(field.name())) {
      Log::error() << "This FieldSet3D variable " << field.name() << " does not match "
                         << "corresponding serialized FieldSet3D variable." << std::endl;
      throw eckit::BadParameter("FieldSet3D::deserialize failed: different variable name", Here());
    }
    const int other_field_shape_0 = vect[index++];
    const int other_field_shape_1 = vect[index++];
    if (other_field_shape_0 != field.shape(0) || other_field_shape_1 != field.shape(1)) {
      Log::error() << "This FieldSet3D variable " << field.name() << " sizes: "
                         << field.shape(0) << ", " << field.shape(1) << std::endl;
      Log::error() << "Serialized FieldSet3D variable " << field.name() << " sizes: "
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

// -----------------------------------------------------------------------------

void FieldSet3D::removeFields(const Variables & vars) {
  util::removeFieldsFromFieldSet(fset_, vars.variables());
}

// -----------------------------------------------------------------------------

bool FieldSet3D::compare_with(const FieldSet3D & other,
                              const double & tol,
                              const util::ToleranceType & tolType) const {
  return util::compareFieldSets(comm_, fset_, other.fset_, tol, tolType);
}

// -----------------------------------------------------------------------------

std::string FieldSet3D::getGridUid() const {
  return util::getGridUid(fset_);
}

// -----------------------------------------------------------------------------

}  // namespace oops
