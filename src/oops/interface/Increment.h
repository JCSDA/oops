/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2020 UCAR.
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

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Serializable.h"
#include "oops/util/Timer.h"

namespace oops {

namespace interface {

/// Increment: Difference between two model states.
/// Some fields that are present in a State may not be present in an Increment.
template <typename MODEL>
class Increment : public oops::GeneralizedDepartures,
                  public util::Serializable,
                  private util::ObjectCounter<Increment<MODEL> > {
  typedef typename MODEL::Increment  Increment_;
  typedef oops::Geometry<MODEL>      Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "oops::Increment";}

  /// Constructor for specified \p geometry, with \p variables, valid on \p date
  Increment(const Geometry_ & geometry, const Variables & variables, const util::DateTime & date);
  /// Copies \p other increment, changing its resolution to \p geometry
  Increment(const Geometry_ & geometry, const Increment & other);
  /// Creates Increment with the same geometry and variables as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero increment
  Increment(const Increment &, const bool copy = true);

  /// Destructor (defined explicitly for timing and tracing)
  virtual ~Increment();

  /// Set this Increment to be difference between \p state1 and \p state2
  void diff(const State_ & state1, const State_ & state2);

  /// Accessor to the time of this Increment
  const util::DateTime validTime() const {return increment_->validTime();}
  /// Updates this Increment's valid time by \p dt (used in PseudoModel)
  void updateTime(const util::Duration & dt) {increment_->updateTime(dt);}

  /// Zero out this Increment
  void zero();
  /// Zero out this Increment and set its date to \p date
  void zero(const util::DateTime & date);
  /// Set this Increment to ones (used in tests)
  void ones();
  /// Set Increment according to the configuration (used in Dirac application)
  void dirac(const eckit::Configuration &);

  /// Assignment operator
  Increment & operator =(const Increment &);
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
  void read(const eckit::Configuration &);
  /// Write this Increment out to file
  void write(const eckit::Configuration &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Get local (at \p iter local volume) increment (used in LocalEnsembleSolver)
  LocalIncrement getLocal(const GeometryIterator_ & iter) const;
  /// Set local (at \p iter local volume) increment to be \p gp (used in LocalEnsembleSolver)
  void setLocal(const LocalIncrement & gp, const GeometryIterator_ & iter);

  /// Accessor to geometry associated with this Increment
  Geometry_ geometry() const;

  /// ATLAS FieldSet (used in SABER)
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

  /// Serialize and deserialize (used in 4DEnVar, weak-constraint 4DVar and Block-Lanczos minimizer)
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  /// Interfacing with other oops classes (returns reference to the implementation object)
  const Increment_ & increment() const {return *this->increment_;}
  Increment_ & increment() {return *this->increment_;}

 protected:
  std::unique_ptr<Increment_> increment_;   /// pointer to the Increment implementation

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
Increment<MODEL>::Increment(const Geometry_ & resol, const Increment & other)
  : increment_()
{
  Log::trace() << "Increment<MODEL>::Increment starting" << std::endl;
  util::Timer timer(classname(), "Increment");
  increment_.reset(new Increment_(resol.geometry(), *other.increment_));
  this->setObjectSize(increment_->serialSize()*sizeof(double));
  Log::trace() << "Increment<MODEL>::Increment done" << std::endl;
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
  Log::trace() << "Increment<MODEL>::~Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::diff(const State_ & x1, const State_ & x2) {
  Log::trace() << "Increment<MODEL>::diff starting" << std::endl;
  util::Timer timer(classname(), "diff");
  increment_->diff(x1.state(), x2.state());
  Log::trace() << "Increment<MODEL>::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero() {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  increment_->zero();
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::zero(const util::DateTime & tt) {
  Log::trace() << "Increment<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  increment_->zero(tt);
  Log::trace() << "Increment<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::ones() {
  Log::trace() << "Increment<MODEL>::ones starting" << std::endl;
  util::Timer timer(classname(), "ones");
  increment_->ones();
  Log::trace() << "Increment<MODEL>::ones done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::dirac(const eckit::Configuration & config) {
  Log::trace() << "Increment<MODEL>::dirac starting" << std::endl;
  util::Timer timer(classname(), "dirac");
  increment_->dirac(config);
  Log::trace() << "Increment<MODEL>::dirac done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *increment_ = *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator+=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator+= starting" << std::endl;
  util::Timer timer(classname(), "operator+=");
  *increment_ += *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator-=(const Increment & rhs) {
  Log::trace() << "Increment<MODEL>::operator-= starting" << std::endl;
  util::Timer timer(classname(), "operator-=");
  *increment_ -= *rhs.increment_;
  Log::trace() << "Increment<MODEL>::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> & Increment<MODEL>::operator*=(const double & zz) {
  Log::trace() << "Increment<MODEL>::operator*= starting" << std::endl;
  util::Timer timer(classname(), "operator*=");
  *increment_ *= zz;
  Log::trace() << "Increment<MODEL>::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::axpy(const double & zz, const Increment & dx, const bool check) {
  Log::trace() << "Increment<MODEL>::axpy starting" << std::endl;
  util::Timer timer(classname(), "axpy");
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
  increment_->schur_product_with(*dx.increment_);
  Log::trace() << "Increment<MODEL>::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::random() {
  Log::trace() << "Increment<MODEL>::random starting" << std::endl;
  util::Timer timer(classname(), "random");
  increment_->random();
  Log::trace() << "Increment<MODEL>::random done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::accumul(const double & zz, const State_ & xx) {
  Log::trace() << "Increment<MODEL>::accumul starting" << std::endl;
  util::Timer timer(classname(), "accumul");
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
  increment_->setLocal(gp, iter.geometryiter());
  Log::trace() << "Increment<MODEL>::setLocal done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "Increment<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  increment_->read(conf);
  Log::trace() << "Increment<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "Increment<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  increment_->write(conf);
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
oops::Geometry<MODEL> Increment<MODEL>::geometry() const {
  Log::trace() << "Increment<MODEL>::geometry starting" << std::endl;
  util::Timer timer(classname(), "geometry");
  oops::Geometry<MODEL> geom(increment_->geometry());
  Log::trace() << "Increment<MODEL>::geometry done" << std::endl;
  return geom;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::setAtlas(atlas::FieldSet * atlasFieldSet) const {
  Log::trace() << "Increment<MODEL>::setAtlas starting" << std::endl;
  util::Timer timer(classname(), "setAtlas");
  increment_->setAtlas(atlasFieldSet);
  Log::trace() << "Increment<MODEL>::setAtlas done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::toAtlas(atlas::FieldSet * atlasFieldSet) const {
  Log::trace() << "Increment<MODEL>::toAtlas starting" << std::endl;
  util::Timer timer(classname(), "toAtlas");
  increment_->toAtlas(atlasFieldSet);
  Log::trace() << "Increment<MODEL>::toAtlas done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::fromAtlas(atlas::FieldSet * atlasFieldSet) {
  Log::trace() << "Increment<MODEL>::fromAtlas starting" << std::endl;
  util::Timer timer(classname(), "fromAtlas");
  increment_->fromAtlas(atlasFieldSet);
  Log::trace() << "Increment<MODEL>::fromAtlas done" << std::endl;
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

}  // namespace interface

// -----------------------------------------------------------------------------
/// \brief Increment class used in oops
///
/// \details
/// Adds extra methods that do not need to be implemented in the implementations:
/// - timeComm()  (accessor to the MPI communicator in time - collection of processes
///                holding the data needed to represent the state in a particular region
///                of space X_i and throughout the whole time interval for which DA is done)
/// - variables() (accessor to variables in this Increment)
/// - shift_forward
/// - shift_backward
/// - toAtlas, atlas
///
/// Adds communication through time to the following Increment methods:
/// - dot_product_with
/// - norm
/// - print

template <typename MODEL>
class Increment : public interface::Increment<MODEL> {
  typedef Geometry<MODEL>            Geometry_;

 public:
  /// Constructor for specified \p geometry, with \p variables, valid on \p date
  Increment(const Geometry_ & geometry, const Variables & variables, const util::DateTime & date);
  /// Copies \p other increment, changing its resolution to \p geometry
  Increment(const Geometry_ & geometry, const Increment & other);
  /// Creates Increment with the same geometry and variables as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero increment
  Increment(const Increment & other, const bool copy = true);

  /// Accessor to the time communicator
  const eckit::mpi::Comm & timeComm() const {return *timeComm_;}
  /// Accessor to Variables stored in this increment
  const Variables & variables() const {return variables_;}

  /// Shift forward in time by \p dt
  void shift_forward(const util::DateTime & dt);
  /// Shift backward in time by \p dt
  void shift_backward(const util::DateTime & dt);

  /// Set ATLAS fieldset associated with this Increment internally
  void toAtlas();
  /// Allow to access base class's method as well
  using interface::Increment<MODEL>::toAtlas;
  /// Accessors to the ATLAS fieldset
  atlas::FieldSet & atlas() {return atlasFieldSet_;}
  const atlas::FieldSet & atlas() const {return atlasFieldSet_;}

  /// dot product with the \p other increment
  double dot_product_with(const Increment & other) const;
  /// Norm for diagnostics
  double norm() const;

 private:
  void print(std::ostream &) const override;

  Variables variables_;                /// Variables stored in this Increment
  const eckit::mpi::Comm * timeComm_;  /// pointer to the MPI communicator in time
  atlas::FieldSet atlasFieldSet_;      /// Atlas fields associated with this Increment
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & geometry, const Variables & variables,
                            const util::DateTime & date):
  interface::Increment<MODEL>(geometry, variables, date), variables_(variables),
  timeComm_(&geometry.timeComm())
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Geometry_ & geometry, const Increment & other):
  interface::Increment<MODEL>(geometry, other), variables_(other.variables_),
  timeComm_(other.timeComm_)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Increment<MODEL>::Increment(const Increment & other, const bool copy):
  interface::Increment<MODEL>(other, copy), variables_(other.variables_),
  timeComm_(other.timeComm_)
{}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::dot_product_with(const Increment & dx) const {
  double zz = interface::Increment<MODEL>::dot_product_with(dx);
  timeComm_->allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double Increment<MODEL>::norm() const {
  double zz = interface::Increment<MODEL>::norm();
  zz *= zz;
  timeComm_->allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  zz = sqrt(zz);
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::shift_forward(const util::DateTime & begin) {
  Log::trace() << "Increment<MODEL>::Increment shift_forward starting" << std::endl;
  static int tag = 159357;
  size_t mytime = timeComm_->rank();

// Send values of M.dx_i at end of my subwindow to next subwindow
  if (mytime + 1 < timeComm_->size()) {
    oops::mpi::send(*timeComm_, *this, mytime+1, tag);
  }

// Receive values at beginning of my subwindow from previous subwindow
  if (mytime > 0) {
    oops::mpi::receive(*timeComm_, *this, mytime-1, tag);
  } else {
    this->zero(begin);
  }

  ++tag;
  Log::trace() << "Increment<MODEL>::Increment shift_forward done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::shift_backward(const util::DateTime & end) {
  Log::trace() << "Increment<MODEL>::Increment shift_backward starting" << std::endl;
  static int tag = 30951;
  size_t mytime = timeComm_->rank();

// Send values of dx_i at start of my subwindow to previous subwindow
  if (mytime > 0) {
    oops::mpi::send(*timeComm_, *this, mytime-1, tag);
  }

// Receive values at end of my subwindow from next subwindow
  if (mytime + 1 < timeComm_->size()) {
    oops::mpi::receive(*timeComm_, *this, mytime+1, tag);
  } else {
    this->zero(end);
  }

  ++tag;
  Log::trace() << "Increment<MODEL>::Increment shift_backward done" << std::endl;
}


// -----------------------------------------------------------------------------
template<typename MODEL>
void Increment<MODEL>::toAtlas() {
  interface::Increment<MODEL>::toAtlas(&atlasFieldSet_);
  this->increment_.reset();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Increment<MODEL>::print(std::ostream & os) const {
  if (timeComm_->size() > 1) {
    gatherPrint(os, this->increment(), *timeComm_);
  } else {
    os << this->increment();
  }
}

// -----------------------------------------------------------------------------
/// Add on \p dx incrment to model state \p xx
template <typename MODEL>
State<MODEL> & operator+=(State<MODEL> & xx, const Increment<MODEL> & dx) {
  Log::trace() << "operator+=(State, Increment) starting" << std::endl;
  util::Timer timer("oops::Increment", "operator+=(State, Increment)");
  xx.state() += dx.increment();
  Log::trace() << "operator+=(State, Increment) done" << std::endl;
  return xx;
}


}  // namespace oops

#endif  // OOPS_INTERFACE_INCREMENT_H_
