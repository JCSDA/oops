/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_OBSERRORBASE_H_
#define OOPS_GENERIC_OBSERRORBASE_H_

#include <map>
#include <memory>
#include <string>

#include <boost/noncopyable.hpp>
#include "eckit/config/Configuration.h"

#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Base class for generic implementations of observation error covariance matrices.
///
/// Use this class as a base class for generic implementations,
/// and interface::ObsErrorBase as a base class for OBS-specific implementations.
///
/// Note: generic implementations need to provide a constructor with the following signature:
///
///     ObsErrorBase(const eckit::Configuration &config, ObsSpace<OBS> &obsspace);
template<typename OBS>
class ObsErrorBase : public util::Printable,
                     private boost::noncopyable {
  typedef ObsVector<OBS>        ObsVector_;
  typedef ObsSpace<OBS>         ObsSpace_;

 public:
  ObsErrorBase() = default;
  virtual ~ObsErrorBase() = default;

/// Multiply a Departure \p dy by \f$R\f$.
  virtual void multiply(ObsVector_ & dy) const = 0;
/// Multiply a Departure \p dy by \f$R^{-1}\f$.
  virtual void inverseMultiply(ObsVector_ & dy) const = 0;

/// Generate random perturbation in \p dy.
  virtual void randomize(ObsVector_ & dy) const = 0;

/// Save obs errors to the group \p name.
  virtual void save(const std::string &name) const = 0;

/// Return a copy of obs error std. dev. If this ObsVector_ is modified (e.g. by obs filters),
/// it should be passed back to update() to ensure the covariance matrix stays consistent.
  virtual ObsVector_ obserrors() const = 0;

/// Set the diagonal of the covariance matrix to \p stddev squared.
  virtual void update(const ObsVector_ &stddev) = 0;

/// Return the vector of inverse obs error variances.
  virtual ObsVector_ inverseVariance() const = 0;

/// Get mean error for Jo table.
  virtual double getRMSE() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

/// A factory creating instances of concrete subclasses of ObsErrorBase.
template <typename OBS>
class ObsErrorFactory {
  typedef ObsErrorBase<OBS> ObsErrorBase_;
  typedef ObsSpace<OBS>     ObsSpace_;
 public:
  static std::unique_ptr<ObsErrorBase_> create(const eckit::Configuration &,
                                                    const ObsSpace_ &);
  virtual ~ObsErrorFactory() = default;
 protected:
  explicit ObsErrorFactory(const std::string &);
 private:
  virtual std::unique_ptr<ObsErrorBase_> make(const eckit::Configuration &, const ObsSpace_ &) = 0;
  static std::map < std::string, ObsErrorFactory<OBS> * > & getMakers() {
    static std::map < std::string, ObsErrorFactory<OBS> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of ObsErrorFactory able to create instances of T (a concrete subclass of
/// ObsErrorBase<OBS>). Passes ObsSpace<OBS> to the constructor of T.
template<class OBS, class T>
class ObsErrorMaker : public ObsErrorFactory<OBS> {
  typedef ObsErrorBase<OBS> ObsErrorBase_;
  typedef ObsSpace<OBS> ObsSpace_;
  std::unique_ptr<ObsErrorBase_> make(const eckit::Configuration & conf,
                                      const ObsSpace_ & obs) override
    { return std::make_unique<T>(conf, obs); }
 public:
  explicit ObsErrorMaker(const std::string & name) : ObsErrorFactory<OBS>(name) {}
};

// =============================================================================

template <typename OBS>
ObsErrorFactory<OBS>::ObsErrorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    throw std::runtime_error(name + " already registered in obs error factory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename OBS>
std::unique_ptr<ObsErrorBase<OBS>>
ObsErrorFactory<OBS>::create(const eckit::Configuration & conf, const ObsSpace_ & obs) {
  Log::trace() << "ObsErrorFactory<OBS>::create starting" << std::endl;
  const std::string id = conf.getString("covariance model", "diagonal");
  Log::trace() << "ObsError matrix type is: " << id << std::endl;
  typename std::map<std::string, ObsErrorFactory<OBS>*>::iterator
    jerr = getMakers().find(id);
  if (jerr == getMakers().end()) {
    throw std::runtime_error(id + " does not exist in obs error factory.");
  }
  std::unique_ptr<ObsErrorBase<OBS>> ptr(jerr->second->make(conf, obs));
  Log::trace() << "ObsErrorFactory<OBS>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_OBSERRORBASE_H_
