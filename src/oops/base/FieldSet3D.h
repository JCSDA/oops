/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iterator>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Wrapper around atlas::FieldSet, storing the valid time associated with the fieldset.
/// FieldSet3D implements all methods that are required to be implemented in the DATA
/// associated with oops::DataSetBase<DATA,...> (base class for the 4D ensemble storage).
/// Note: FieldSet3D assumes that the valid time of the data does not change during
/// the object's lifecycle, unlike oops::State and oops::Increment.
/// FieldSet3D assumes that the variables of the data can change during its lifecycle.
class FieldSet3D : public util::Serializable,
                   public util::Printable {
 public:
  using iterator       = std::vector<atlas::Field>::iterator;
  using const_iterator = std::vector<atlas::Field>::const_iterator;

  /// Constructors
  /// @brief Creates FieldSet3D with an empty atlas::FieldSet
  FieldSet3D(const util::DateTime &, const eckit::mpi::Comm &);
  /// @brief Copy constructor (deep copy)
  FieldSet3D(const FieldSet3D &);

  /// Initialization of atlas::FieldSet
  /// @brief Initialize atlas::FieldSet only, consistency check for valid time
  void allocateOnly(const atlas::FieldSet &);
  void allocateOnly(const FieldSet3D & other);
  /// @brief Initialize atlas::FieldSet with a deep copy, consistency check for valid time
  void deepCopy(const atlas::FieldSet &);
  void deepCopy(const FieldSet3D & other);
  /// @brief Initialize atlas::FieldSet with a shallow copy, consistency check for valid time
  void shallowCopy(const atlas::FieldSet &);
  void shallowCopy(const FieldSet3D & other);
  /// @brief Initialize atlas::FieldSet, possibly with a constant value
  void init(const atlas::FunctionSpace &,
            const oops::Variables &);
  void init(const atlas::FunctionSpace &,
            const oops::Variables &,
            const double & value);
  /// @brief Initialize atlas::FieldSet with random fields
  void randomInit(const atlas::FunctionSpace &,
                  const oops::Variables &);

  /// Accessors
  /// @brief Return variables
  const Variables & variables() const;
  /// @brief Return valid time
  const util::DateTime validTime() const {return validTime_;}
  util::DateTime & validTime() {return validTime_;}
  /// @brief Return communicator
  const eckit::mpi::Comm & commGeom() const {return comm_;}
  /// @brief Return atlas::FieldSet
  const atlas::FieldSet & fieldSet() const {return fset_;}
  atlas::FieldSet & fieldSet() {return fset_;}
  /// @brief Return name
  const std::string & name() const {return name_ == "" ? fset_.name() : name_;}
  std::string & name() {return name_;}
  const atlas::Field & operator[](const int & fieldIndex) const {return fset_[fieldIndex];}
  atlas::Field & operator[](const int & fieldIndex) {return fset_[fieldIndex];}
  const atlas::Field & operator[](const std::string & fieldName) const {return fset_[fieldName];}
  atlas::Field & operator[](const std::string & fieldName) {return fset_[fieldName];}
  const atlas::Field & operator[](const Variable & var) const {return fset_[var.name()];}
  atlas::Field & operator[](const Variable & var) {return fset_[var.name()];}


  /// Arithmetic operations
  void zero();
  FieldSet3D & operator+=(const FieldSet3D &);
  FieldSet3D & operator-=(const FieldSet3D &);
  FieldSet3D & operator*=(const FieldSet3D &);
  FieldSet3D & operator*=(const double &);
  FieldSet3D & operator/=(const FieldSet3D &);
  /// @brief Computes dot product of this FieldSet3D with the \p other FieldSet3D
  ///        only for specified variables \p vars.
  double dot_product_with(const FieldSet3D &, const Variables &) const;
  double norm(const Variables &) const;
  void sqrt();

  /// Read / write
  void read(const atlas::FunctionSpace &,
            const oops::Variables &,
            const eckit::LocalConfiguration &);
  void write(const eckit::LocalConfiguration &) const;

  /// Serialize / deserializee
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;


  /// Utilities
  void removeFields(const Variables & vars);
  bool compare_with(const FieldSet3D &,
                    const double & tol = 1.0e-12,
                    const util::ToleranceType & tolType = util::ToleranceType::absolute) const;
  std::vector<std::string> field_names() const {return fset_.field_names();}
  size_t size() const {return fset_.size();}
  bool empty() const {return fset_.empty();}
  std::string getGridUid() const;
  bool has(const std::string & fieldName) const {return fset_.has(fieldName);}
  void add(const atlas::Field & field) {fset_.add(field);}
  void clear() {fset_.clear();}
  iterator begin() {return fset_.begin();}
  iterator end() {return fset_.end();}
  const_iterator begin() const {return fset_.begin();}
  const_iterator end() const {return fset_.end();}
  const_iterator cbegin() const {return fset_.cbegin();}
  const_iterator cend() const {return fset_.cend();}
  const atlas::field::FieldSetImpl * get() const {return fset_.get();}

 private:
  oops::Variables currentVariables() const;
  void print(std::ostream &) const override;

  atlas::FieldSet fset_;
  util::DateTime validTime_;
  const eckit::mpi::Comm & comm_;
  std::string name_;
  mutable oops::Variables vars_;
};

// -----------------------------------------------------------------------------

/// @brief Initializes FieldSet3D to have the same fields as the \p other; values
///        are not copied.
FieldSet3D initFieldSet3D(const FieldSet3D &);
FieldSet3D initFieldSet3D(const util::DateTime &,
                          const eckit::mpi::Comm &,
                          const atlas::FieldSet &);

/// @brief Initializes FieldSet3D to have the same fields as the \p other; values
///        are copied, through a deep copy (default) or a shallow copy
FieldSet3D copyFieldSet3D(const FieldSet3D &,
                          const bool & shallow = false);
FieldSet3D copyFieldSet3D(const util::DateTime &,
                          const eckit::mpi::Comm &,
                          const atlas::FieldSet &,
                          const bool & shallow = false);

/// @brief Initializes FieldSet3D with random fields
FieldSet3D randomFieldSet3D(const util::DateTime &,
                            const eckit::mpi::Comm &,
                            const atlas::FunctionSpace &,
                            const oops::Variables &);
}  // namespace oops

