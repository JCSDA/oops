/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_MODELAUXCONTROL_H_
#define OOPS_INTERFACE_MODELAUXCONTROL_H_

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Auxiliary state related to model (could be e.g. model bias), not used at the moment.
/// \details
/// This class is used to manipulate parameters of the model that can be estimated in the
/// assimilation. This includes model bias but could be used for other parameters.
/// This is sometimes referred to as augmented state or augmented control variable in the
/// literature.
/// This class calls the model's implementation of ModelAuxControl.

// -----------------------------------------------------------------------------

template <typename MODEL>
class ModelAuxControl : public util::Printable,
                        public util::Serializable,
                        private util::ObjectCounter<ModelAuxControl<MODEL> > {
  typedef typename MODEL::ModelAuxControl      ModelAuxControl_;
  typedef Geometry<MODEL>            Geometry_;

 public:
  static const std::string classname() {return "oops::ModelAuxControl";}

  ModelAuxControl(const Geometry_ &, const eckit::Configuration &);
  /// Copies \p other ModelAuxControl, changing its resolution to \p resol
  ModelAuxControl(const Geometry_ & resol, const ModelAuxControl & other);
  /// Creates ModelAuxControl with the same structure as \p other.
  /// Copies \p other if \p copy is true, otherwise creates zero ModelAuxControl
  explicit ModelAuxControl(const ModelAuxControl &, const bool copy = true);
  /// Destructor (defined explicitly for timing and tracing)
  ~ModelAuxControl();

  /// const Accessor
  const ModelAuxControl_ & modelauxcontrol() const {return *aux_;}
  /// Accessor
  ModelAuxControl_ & modelauxcontrol() {return *aux_;}

  /// Read this ModelAuxControl from file
  void read(const eckit::Configuration &);
  /// Write this ModelAuxControl out to file
  void write(const eckit::Configuration &) const;
  /// Norm (used in tests)
  double norm() const;

  /// Serialize and deserialize
  size_t serialSize() const override {return 0;}
  void serialize(std::vector<double> &) const override {}
  void deserialize(const std::vector<double> &, size_t &) override {}

  /// Assignment
  ModelAuxControl & operator=(const ModelAuxControl &);

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<ModelAuxControl_> aux_;
};

// =============================================================================

template<typename MODEL>
ModelAuxControl<MODEL>::ModelAuxControl(const Geometry_ & resol,
                                        const eckit::Configuration & conf) : aux_()
{
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxControl");
  aux_.reset(new ModelAuxControl_(resol.geometry(), conf));
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxControl<MODEL>::ModelAuxControl(const Geometry_ & resol,
                                        const ModelAuxControl & other) : aux_()
{
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl interpolated starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxControl");
  aux_.reset(new ModelAuxControl_(resol.geometry(), *other.aux_));
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxControl<MODEL>::ModelAuxControl(const ModelAuxControl & other,
                                        const bool copy) : aux_()
{
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl copy starting" << std::endl;
  util::Timer timer(classname(), "ModelAuxControl");
  aux_.reset(new ModelAuxControl_(*other.aux_, copy));
  Log::trace() << "ModelAuxControl<MODEL>::ModelAuxControl copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxControl<MODEL>::~ModelAuxControl() {
  Log::trace() << "ModelAuxControl<MODEL>::~ModelAuxControl starting" << std::endl;
  util::Timer timer(classname(), "~ModelAuxControl");
  aux_.reset();
  Log::trace() << "ModelAuxControl<MODEL>::~ModelAuxControl done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxControl<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "ModelAuxControl<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  aux_->read(conf);
  Log::trace() << "ModelAuxControl<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxControl<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ModelAuxControl<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  aux_->write(conf);
  Log::trace() << "ModelAuxControl<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double ModelAuxControl<MODEL>::norm() const {
  Log::trace() << "ModelAuxControl<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = aux_->norm();
  Log::trace() << "ModelAuxControl<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ModelAuxControl<MODEL> & ModelAuxControl<MODEL>::operator=(const ModelAuxControl & rhs) {
  Log::trace() << "ModelAuxControl<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *aux_ = *rhs.aux_;
  Log::trace() << "ModelAuxControl<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ModelAuxControl<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ModelAuxControl<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *aux_;
  Log::trace() << "ModelAuxControl<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_MODELAUXCONTROL_H_
