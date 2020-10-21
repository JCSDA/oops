/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <boost/noncopyable.hpp>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

/*! \brief Base class for Generic interpolation
 *
 * \details the **oops::InterpolatorBase** class provides a generic framework
 * for interpolation.
 *
 */

class InterpolatorBase : public util::Printable,
                         private boost::noncopyable {
 public:
  virtual ~InterpolatorBase() {}

  virtual void apply(const atlas::FieldSet &, atlas::FieldSet &) = 0;
  virtual void apply(const atlas::Field &, atlas::Field &) = 0;

  virtual void apply_ad(const atlas::FieldSet &, atlas::FieldSet &) = 0;

  virtual int write(const eckit::Configuration &) {return 1;}

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------
/// InterpolatorFactory: Factory for creating Interpolator objects

class InterpolatorFactory {
 public:
  static InterpolatorBase * create(const eckit::Configuration &,
                              const atlas::FunctionSpace &,
                              const atlas::FunctionSpace &,
                              const atlas::field::FieldSetImpl * = nullptr);
  virtual ~InterpolatorFactory() = default;
 protected:
  explicit InterpolatorFactory(const std::string &);
 private:
  virtual InterpolatorBase * make(const eckit::Configuration &,
                                  const atlas::FunctionSpace &,
                                  const atlas::FunctionSpace &,
                                  const atlas::field::FieldSetImpl *) = 0;
  static std::map < std::string, InterpolatorFactory * > & getMakers() {
    static std::map < std::string, InterpolatorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------
/// Maker class for Interpolator Factory
///
/// Used to instantiate concrete Interpolator objects of
/// derived type T

template<class T>
class InterpolatorMaker : public InterpolatorFactory {
  virtual InterpolatorBase * make(const eckit::Configuration & conf,
                                  const atlas::FunctionSpace & fs1,
                                  const atlas::FunctionSpace & fs2,
                                  const atlas::field::FieldSetImpl * masks)
    { return new T(conf, fs1, fs2, masks); }
 public:
  explicit InterpolatorMaker(const std::string & name)
    : InterpolatorFactory(name) {}
};



// -----------------------------------------------------------------------------

}  // namespace oops
