/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_MODELAUXINCREMENT_H_
#define OOPS_INTERFACE_MODELAUXINCREMENT_H_

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Auxiliary Increment related to model, not used at the moment.
/// \details
/// This class calls the model's implementation of ModelAuxIncrement.
// -----------------------------------------------------------------------------

/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of Parameters.
///
/// In the former case, they should provide constructors with the following signatures:
///
///    ModelAuxIncrement(const Geometry_ &, const eckit::Configuration &);
///    ModelAuxIncrement(const ModelAuxIncrement &, const eckit::Configuration &);
///
/// In the latter case, the implementer should first define a subclass of Parameters holding the
/// settings of the ModelAuxIncrement implementation in question. That implementation should then
/// typedef `Parameters_` to the name of that subclass and provide constructors with the following
/// signatures:
///
///    ModelAuxIncrement(const Geometry_ &, const Parameters_ &);
///    ModelAuxIncrement(const ModelAuxIncrement &, const Parameters_ &);
///
/// The implementations of the ModelAuxIncrement and ModelAuxCovariance interfaces for the same
/// MODEL should use the same Parameters subclass.
template <typename MODEL>
class ModelAuxIncrement : public util::Printable,
                          public util::Serializable,
                          private util::ObjectCounter<ModelAuxIncrement<MODEL> > {
  typedef typename MODEL::ModelAuxIncrement     ModelAuxIncrement_;
  typedef Geometry<MODEL>             Geometry_;
  typedef ModelAuxControl<MODEL>      ModelAuxControl_;

 public:
  /// Set to ModelAuxIncrement_::Parameters_ if ModelAuxIncrement_ provides a type called
  /// Parameters_ and to GenericParameters (a thin wrapper of an eckit::LocalConfiguration object)
  /// if not.
  typedef TParameters_IfAvailableElseFallbackType_t<ModelAuxIncrement_, GenericParameters>
    Parameters_;

  static const std::string classname() {return "oops::ModelAuxIncrement";}

  /// Constructor for specified \p resol and \p conf
  ModelAuxIncrement(const Geometry_ & resol, const eckit::Configuration & conf);
  ModelAuxIncrement(const Geometry_ &, const Parameters_ &);
  /// Copies \p other ModelAuxIncrement if \p copy is true,
  /// otherwise creates zero ModelAuxIncrement with same variables and geometry
  explicit ModelAuxIncrement(const ModelAuxIncrement & other, const bool copy = true);
  /// Copies \p other ModelAuxIncrement, reading extra information from \p conf
  ModelAuxIncrement(const ModelAuxIncrement & other, const eckit::Configuration & conf);
  ModelAuxIncrement(const ModelAuxIncrement &, const Parameters_ &);
  /// Destructor (defined explicitly for timing and tracing)
  ~ModelAuxIncrement();

  /// const Accessor
  const ModelAuxIncrement_ & modelauxincrement() const {return *aux_;}
  /// Accessor
  ModelAuxIncrement_ & modelauxincrement() {return *aux_;}

  /// Sets this ModelAuxIncrement to the difference between two ModelAuxControl objects
  void diff(const ModelAuxControl_ &, const ModelAuxControl_ &);
  /// Zero out this ModelAuxIncrement
  void zero();
  /// Linear algebra operators
  ModelAuxIncrement & operator=(const ModelAuxIncrement &);
  ModelAuxIncrement & operator+=(const ModelAuxIncrement &);
  ModelAuxIncrement & operator-=(const ModelAuxIncrement &);
  ModelAuxIncrement & operator*=(const double &);
  void axpy(const double &, const ModelAuxIncrement &);
  /// dot product with the \p other ModelAuxIncrement
  double dot_product_with(const ModelAuxIncrement & other) const;

  /// Read this ModelAuxIncrement from file
  void read(const eckit::Configuration &);
  /// Write this ModelAuxIncrement out to file
  void write(const eckit::Configuration &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Serialize and deserialize (used in 4DEnVar, weak-constraint 4DVar and Block-Lanczos minimizer)
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ModelAuxIncrement_> aux_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ModelAuxControl<MODEL> & operator+=(ModelAuxControl<MODEL> & xx,
                                    const ModelAuxIncrement<MODEL> & dx) {
  Log::trace() << "operator+=(ModelAuxControl, ModelAuxIncrement) starting" << std::endl;
  util::Timer timer("oops::ModelAuxIncrement", "operator+=ModelAuxControl");
  xx.modelauxcontrol() += dx.modelauxincrement();
  Log::trace() << "operator+=(ModelAuxControl, ModelAuxIncrement) done" << std::endl;
  return xx;
}

// =============================================================================

template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const Geometry_ & resol,
                                            const Parameters_ & parameters) : aux_()
{
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrement");
  aux_.reset(new ModelAuxIncrement_(
               resol.geometry(),
               parametersOrConfiguration<HasParameters_<ModelAuxIncrement_>::value>(parameters)));
  this->setObjectSize(aux_->serialSize()*sizeof(double));
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const Geometry_ & resol,
                                            const eckit::Configuration & conf)
  : ModelAuxIncrement(resol, validateAndDeserialize<Parameters_>(conf))
{}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const ModelAuxIncrement & other,
                                            const bool copy) : aux_()
{
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement copy starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrement");
  aux_.reset(new ModelAuxIncrement_(*other.aux_, copy));
  this->setObjectSize(aux_->serialSize()*sizeof(double));
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement copy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const ModelAuxIncrement & other,
                                            const Parameters_ & parameters) : aux_()
{
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement interpolated starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxIncrement");
  aux_.reset(new ModelAuxIncrement_(
               *other.aux_,
               parametersOrConfiguration<HasParameters_<ModelAuxIncrement_>::value>(parameters)));
  this->setObjectSize(aux_->serialSize()*sizeof(double));
  Log::trace() << "ModelAuxIncrement<MODEL>::ModelAuxIncrement interpolated done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::ModelAuxIncrement(const ModelAuxIncrement & other,
                                            const eckit::Configuration & conf)
  : ModelAuxIncrement(other, validateAndDeserialize<Parameters_>(conf))
{}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL>::~ModelAuxIncrement() {
  Log::trace() << "ModelAuxIncrement<MODEL>::~ModelAuxIncrement starting" << std::endl;
  util::Timer timer(classname(), "~ModelAuxIncrement");
  aux_.reset();
  Log::trace() << "ModelAuxIncrement<MODEL>::~ModelAuxIncrement done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::diff(const ModelAuxControl_ & x1, const ModelAuxControl_ & x2) {
  Log::trace() << "ModelAuxIncrement<MODEL>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  aux_->diff(x1.modelauxcontrol(), x2.modelauxcontrol());
  Log::trace() << "ModelAuxIncrement<MODEL>::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::zero() {
  Log::trace() << "ModelAuxIncrement<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  aux_->zero();
  Log::trace() << "ModelAuxIncrement<MODEL>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator=(const ModelAuxIncrement & rhs) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *aux_ = *rhs.aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator+=(const ModelAuxIncrement & rhs) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *aux_ += *rhs.aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator-=(const ModelAuxIncrement & rhs) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *aux_ -= *rhs.aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
ModelAuxIncrement<MODEL> & ModelAuxIncrement<MODEL>::operator*=(const double & zz) {
  Log::trace() << "ModelAuxIncrement<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *aux_ *= zz;
  Log::trace() << "ModelAuxIncrement<MODEL>::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::axpy(const double & zz, const ModelAuxIncrement & dx) {
  Log::trace() << "ModelAuxIncrement<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  aux_->axpy(zz, *dx.aux_);
  Log::trace() << "ModelAuxIncrement<MODEL>::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ModelAuxIncrement<MODEL>::dot_product_with(const ModelAuxIncrement & dx) const {
  Log::trace() << "ModelAuxIncrement<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = aux_->dot_product_with(*dx.aux_);
  Log::trace() << "ModelAuxIncrement<MODEL>::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ModelAuxIncrement<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(conf);
  Log::trace() << "ModelAuxIncrement<MODEL>::read done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ModelAuxIncrement<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(conf);
  Log::trace() << "ModelAuxIncrement<MODEL>::write done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
double ModelAuxIncrement<MODEL>::norm() const {
  Log::trace() << "ModelAuxIncrement<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ModelAuxIncrement<MODEL>::norm done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
size_t ModelAuxIncrement<MODEL>::serialSize() const {
  Log::trace() << "ModelAuxIncrement<MODEL>::serialSize" << std::endl;
  util::Timer timer(classname(), "serialSize");
  return aux_->serialSize();
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::serialize(std::vector<double> & vect) const {
  Log::trace() << "ModelAuxIncrement<MODEL>::serialize starting" << std::endl;
  util::Timer timer(classname(), "serialize");
  aux_->serialize(vect);
  Log::trace() << "ModelAuxIncrement<MODEL>::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::deserialize(const std::vector<double> & vect, size_t & current) {
  Log::trace() << "ModelAuxIncrement<MODEL>::deserialize starting" << std::endl;
  util::Timer timer(classname(), "deserialize");
  aux_->deserialize(vect, current);
  Log::trace() << "ModelAuxIncrement<MODEL>::deserialize done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void ModelAuxIncrement<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ModelAuxIncrement<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ModelAuxIncrement<MODEL>::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELAUXINCREMENT_H_
