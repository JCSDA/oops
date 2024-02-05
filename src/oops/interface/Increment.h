/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_INCREMENT_H_
#define OOPS_INTERFACE_INCREMENT_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/Geometry.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/base/WriteParametersBase.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/parameters/HasDiracParameters_.h"
#include "oops/util/parameters/HasReadParameters_.h"
#include "oops/util/parameters/HasWriteParameters_.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/Serializable.h"
#include "oops/util/Timer.h"

namespace oops {

namespace interface {

/// Increment: Difference between two model states.
/// Some fields that are present in a State may not be present in an Increment.
///
/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of Parameters.
///
/// In the former case, they should provide dirac(), read() and write() methods with the
/// following signatures:
///
///     void dirac(const eckit::Configuration &);
///     void read(const eckit::Configuration &);
///     void write(const eckit::Configuration &) const;
///
/// In the latter case, the implementer should first define two subclasses of Parameters, holding
/// the settings needed by the dirac() and read() methods, and a subclass of WriteParametersBase,
/// holding the settings needed by the write() method.
/// The implementation of the Increment interface should then typedef `DiracParameters_`,
/// `ReadParameters_` and `WriteParameters_` to the names of these subclasses and provide dirac(),
/// read() and write() methods with the following signatures:
///
///     void dirac(const DiracParameters_ &);
///     void read(const ReadParameters_ &);
///     void write(const WriteParameters_ &) const;

template <typename MODEL>
class Increment : public oops::GeneralizedDepartures,
                  public util::Serializable,
                  private util::ObjectCounter<Increment<MODEL> > {
  typedef typename MODEL::Increment  Increment_;
  typedef oops::Geometry<MODEL>      Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef oops::State<MODEL>         State_;

 public:
  /// Set to Increment_::DiracParameters_ if Increment_ provides a type called DiracParameters_ and
  /// to GenericParameters (a thin wrapper of an eckit::LocalConfiguration object) if not.
  typedef TDiracParameters_IfAvailableElseFallbackType_t<Increment_, GenericParameters>
    DiracParameters_;
  /// Set to Increment_::ReadParameters_ if Increment_ provides a type called ReadParameters_ and
  /// to GenericParameters (a thin wrapper of an eckit::LocalConfiguration object) if not.
  typedef TReadParameters_IfAvailableElseFallbackType_t<Increment_, GenericParameters>
    ReadParameters_;
  /// Set to Increment_::WriteParameters_ if Increment_ provides a type called WriteParameters_ and
  /// to GenericParameters (a thin wrapper of an eckit::LocalConfiguration object) if not.
  typedef TWriteParameters_IfAvailableElseFallbackType_t<Increment_, GenericWriteParameters>
    WriteParameters_;

  static const std::string classname() {return "oops::Increment";}

 protected:
  /// Constructor for specified \p geometry, with \p variables, valid on \p date
  Increment(const Geometry_ & geometry, const Variables & variables, const util::DateTime & date);
  /// Copies \p other increment, changing its resolution to \p geometry
  /// Uses adjoint resolution change if \p ad is true, though still to the target \p geometry
  Increment(const Geometry_ & geometry, const Increment & other, const bool ad);
  /// Creates Increment with the same geometry and variables as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero increment
  Increment(const Increment &, const bool copy = true);

 public:
  /// Destructor (defined explicitly for timing and tracing)
  virtual ~Increment();

  /// Set this Increment to be difference between \p state1 and \p state2
  void diff(const State_ & state1, const State_ & state2);

  /// Accessor to the time of this Increment
  const util::DateTime validTime() const {return increment_->validTime();}
  /// Updates this Increment's valid time by \p dt (used in PseudoModel)
  void updateTime(const util::Duration & dt) {increment_->updateTime(dt);}
  /// Accessor to variables associated with this Increment
  const Variables & variables() const {return increment_->variables();}

  /// Zero out this Increment
  void zero();
  /// Zero out this Increment and set its date to \p date
  void zero(const util::DateTime & date);
  /// Set this Increment to ones (used in tests)
  void ones();
  /// Set Increment according to the configuration (used in Dirac application)
  void dirac(const DiracParameters_ &);
  /// Set Increment according to the configuration (used in Dirac application)
  void dirac(const eckit::Configuration &);

 protected:
  /// Assignment operator
  Increment & operator =(const Increment &);
 public:
  /// Linear algebra operators
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  /// Add \p w * \p dx to the Increment. If \p check is set, check whether this and \p dx's
  /// dates are the same
  void axpy(const double & w, const Increment & dx, const bool check = true);
  /// Compute dot product of this Increment with \p other
  double dot_product_with(const Increment & other) const;
  /// Compute Schur product of this Increment with \p other, assign to this Increment
  void schur_product_with(const Increment & other);

  /// Randomize the Increment (used in tests)
  void random();
  /// Accumulate (add \p w * \p x to the increment), used in WeightedDiff with Accumulator
  void accumul(const double & w, const State_ & x);

  /// Read this Increment from file
  void read(const ReadParameters_ &);
  /// Read this Increment from file
  void read(const eckit::Configuration &);
  /// Write this Increment out to file
  void write(const WriteParameters_ &) const;
  /// Write this Increment out to file
  void write(const eckit::Configuration &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Get local (at \p iter local volume) increment (used in LocalEnsembleSolver)
  LocalIncrement getLocal(const GeometryIterator_ & iter) const;
  /// Set local (at \p iter local volume) increment to be \p gp (used in LocalEnsembleSolver)
  void setLocal(const LocalIncrement & gp, const GeometryIterator_ & iter);

  /// ATLAS FieldSet interface
  /// For models that are not using ATLAS fieldsets for their own Increment data:
  /// - "toFieldSet" allocates the ATLAS fieldset based on the variables present in the Increment
  ///   and copies Increment data into the fieldset, including halo.
  /// - "fromFieldSet" copies fieldset data back into the Increment (interior points only).
  /// For models that are using ATLAS fieldsets for their own Increment data, fields are shared from
  /// a fieldset to another. A working example is available with the QUENCH testbed of SABER.
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

  /// Serialize and deserialize (used in 4DEnVar, weak-constraint 4DVar and Block-Lanczos minimizer)
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  /// Accessor to MODEL::Increment, used in the other interface classes in oops.
  /// Does not need to be implemented.
  const Increment_ & increment() const {return *this->increment_;}
  Increment_ & increment() {if (fset_) {fset_->clear();} return *this->increment_;}

 protected:
  std::unique_ptr<Increment_> increment_;   /// pointer to the Increment implementation
  mutable std::unique_ptr<FieldSet3D> fset_;

 private:
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & resol, const Variables & vars,
                            const util::DateTime & time)
  : increment_()
{
  Log::trace() << "Increment<MODEL>::Increment starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(resol.geometry(), vars, time));
  this->setObjectSize(increment_->serialSize()*sizeof(double));
  Log::trace() << "Increment<MODEL>::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & resol, const Increment & other, const bool ad)
  : increment_()
{
  Log::trace() << "Increment<MODEL>::Increment chres starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  // Adjoint resolution change currently only used in multiresolution HybridLinearModel.
  // If MODEL attempts to use this without having implemented adjoint resolution change, abort.
  if constexpr (std::is_constructible_v<Increment_,
                const typename MODEL::Geometry &, const Increment_ &, const bool>) {
    increment_.reset(new Increment_(resol.geometry(), *other.increment_, ad));
  } else {
    if (ad) {
      ABORT("Increment<MODEL>::Increment adjoint resolution change requested but not implemented");
    }
    increment_.reset(new Increment_(resol.geometry(), *other.increment_));
  }
  this->setObjectSize(increment_->serialSize()*sizeof(double));
  Log::trace() << "Increment<MODEL>::Increment chres done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::Increment(const Increment & other, const bool copy)
  : increment_()
{
  Log::trace() << "Increment<MODEL>::Increment copy starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(*other.increment_, copy));
  this->setObjectSize(increment_->serialSize()*sizeof(double));
  Log::trace() << "Increment<MODEL>::Increment copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL>::~Increment() {
  Log::trace() << "Increment<MODEL>::~Increment starting" << std::endl;
  util::Timer timer(classname(), "~Increment");
  increment_.reset();
  fset_.reset();
  Log::trace() << "Increment<MODEL>::~Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::diff(const State_ & x1, const State_ & x2) {
  Log::trace() << "Increment<MODEL>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  if (fset_) fset_->clear();
  increment_->diff(x1.state(), x2.state());
  Log::trace() << "Increment<MODEL>::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero() {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  if (fset_) fset_->clear();
  increment_->zero();
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero(const util::DateTime & tt) {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  if (fset_) fset_->clear();
  increment_->zero(tt);
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::ones() {
  Log::trace() << "Increment<MODEL>::ones starting" << std::endl;
  util::Timer timer(classname(), "ones");
  if (fset_) fset_->clear();
  increment_->ones();
  Log::trace() << "Increment<MODEL>::ones done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::dirac(const DiracParameters_ & parameters) {
  dirac(parameters.toConfiguration());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::dirac(const eckit::Configuration & config) {
  Log::trace() << "Increment<MODEL>::dirac starting" << std::endl;
  util::Timer timer(classname(), "dirac");
  if (fset_) fset_->clear();
  increment_->dirac(config);
  Log::trace() << "Increment<MODEL>::dirac done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  if (fset_) fset_->clear();
  *increment_ = *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator+=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  if (fset_) fset_->clear();
  *increment_ += *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator-=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  if (fset_) fset_->clear();
  *increment_ -= *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator*=(const double & zz) {
  Log::trace() << "Increment<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  if (fset_) fset_->clear();
  *increment_ *= zz;
  Log::trace() << "Increment<MODEL>::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::axpy(const double & zz, const Increment & dx, const bool check) {
  Log::trace() << "Increment<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
  if (fset_) fset_->clear();
  increment_->axpy(zz, *dx.increment_, check);
  Log::trace() << "Increment<MODEL>::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::dot_product_with(const Increment & dx) const {
  Log::trace() << "Increment<MODEL>::dot_product_with starting" << std::endl;
  util::Timer timer(classname(), "dot_product_with");
  double zz = increment_->dot_product_with(*dx.increment_);
  Log::trace() << "Increment<MODEL>::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::schur_product_with(const Increment & dx) {
  Log::trace() << "Increment<MODEL>::schur_product_with starting" << std::endl;
  util::Timer timer(classname(), "schur_product_with");
  if (fset_) fset_->clear();
  increment_->schur_product_with(*dx.increment_);
  Log::trace() << "Increment<MODEL>::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::random() {
  Log::trace() << "Increment<MODEL>::random starting" << std::endl;
  util::Timer timer(classname(), "random");
  if (fset_) fset_->clear();
  increment_->random();
  Log::trace() << "Increment<MODEL>::random done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::accumul(const double & zz, const State_ & xx) {
  Log::trace() << "Increment<MODEL>::accumul starting" << std::endl;
  util::Timer timer(classname(), "accumul");
  if (fset_) fset_->clear();
  increment_->accumul(zz, xx.state());
  Log::trace() << "Increment<MODEL>::accumul done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalIncrement Increment<MODEL>::getLocal(const GeometryIterator_ & iter) const {
  Log::trace() << "Increment<MODEL>::getLocal starting" << std::endl;
  util::Timer timer(classname(), "getLocal");
  LocalIncrement gp = increment_->getLocal(iter.geometryiter());
  Log::trace() << "Increment<MODEL>::getLocal done" << std::endl;
  return gp;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment<MODEL>::setLocal(const LocalIncrement & gp,
                                const GeometryIterator_ & iter) {
  Log::trace() << "Increment<MODEL>::setLocal starting" << std::endl;
  util::Timer timer(classname(), "setLocal");
  if (fset_) fset_->clear();
  increment_->setLocal(gp, iter.geometryiter());
  Log::trace() << "Increment<MODEL>::setLocal done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::read(const ReadParameters_ & parameters) {
  increment_->read(parameters.toConfiguration());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::read(const eckit::Configuration & config) {
  Log::trace() << "Increment<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  if (fset_) fset_->clear();
  increment_->read(config);
  Log::trace() << "Increment<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::write(const WriteParameters_ & parameters) const {
  increment_->write(parameters.toConfiguration());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::write(const eckit::Configuration & config) const {
  Log::trace() << "Increment<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  increment_->write(config);
  Log::trace() << "Increment<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::norm() const {
  Log::trace() << "Increment<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = increment_->norm();
  Log::trace() << "Increment<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::toFieldSet(atlas::FieldSet & fset) const {
  Log::trace() << "Increment<MODEL>::toFieldSet starting" << std::endl;
  util::Timer timer(classname(), "toFieldSet");
  increment_->toFieldSet(fset);
  Log::trace() << "Increment<MODEL>::toFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::fromFieldSet(const atlas::FieldSet & fset) {
  Log::trace() << "Increment<MODEL>::fromFieldSet starting" << std::endl;
  util::Timer timer(classname(), "fromFieldSet");
  increment_->fromFieldSet(fset);
  if (fset_) fset_->clear();
  Log::trace() << "Increment<MODEL>::fromFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
size_t Increment<MODEL>::serialSize() const {
  Log::trace() << "Increment<MODEL>::serialSize" << std::endl;
  util::Timer timer(classname(), "serialSize");
  return increment_->serialSize();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::serialize(std::vector<double> & vect) const {
  Log::trace() << "Increment<MODEL>::serialize starting" << std::endl;
  util::Timer timer(classname(), "serialize");
  increment_->serialize(vect);
  Log::trace() << "Increment<MODEL>::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::deserialize(const std::vector<double> & vect, size_t & current) {
  Log::trace() << "Increment<MODEL>::Increment deserialize starting" << std::endl;
  util::Timer timer(classname(), "deserialize");
  if (fset_) fset_->clear();
  increment_->deserialize(vect, current);
  Log::trace() << "Increment<MODEL>::Increment deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Increment<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *increment_;
  Log::trace() << "Increment<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_INCREMENT_H_
