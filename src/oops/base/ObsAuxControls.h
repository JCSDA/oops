/*
 * (C) Copyright 2017-2019 UCAR
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSAUXCONTROLS_H_
#define OOPS_BASE_OBSAUXCONTROLS_H_

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Holds a vector of ObsAuxControl
// -----------------------------------------------------------------------------

template <typename OBS>
class ObsAuxControls : public util::Printable,
                       public util::Serializable {
  typedef ObsAuxControl<OBS>       ObsAuxControl_;
  typedef ObsSpaces<OBS>           ObsSpaces_;
  typedef typename ObsAuxControl_::Parameters_ Parameters_;

 public:
  static const std::string classname() {return "oops::ObsAuxControls";}

  ObsAuxControls(const ObsSpaces_ &, const eckit::Configuration &);
  ObsAuxControls(const ObsSpaces_ &, const std::vector<Parameters_> &);
  explicit ObsAuxControls(const ObsAuxControls &, const bool copy = true);
  ~ObsAuxControls();

/// Access
  std::size_t size() const {return auxs_.size();}
  const ObsAuxControl_ & operator[](const std::size_t ii) const {return *auxs_.at(ii);}
  ObsAuxControl_ & operator[](const std::size_t ii) {return *auxs_.at(ii);}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

  ObsAuxControls & operator=(const ObsAuxControls &);

/// Serialize-Deserialize
  size_t serialSize() const override {return 0;}
  void serialize(std::vector<double> &) const override {}
  void deserialize(const std::vector<double> &, size_t &) override {}

 private:
  void print(std::ostream &) const override;
  std::vector<std::unique_ptr<ObsAuxControl_> > auxs_;
};

// =============================================================================

template<typename OBS>
ObsAuxControls<OBS>::ObsAuxControls(const ObsSpaces_ & odb,
                                    const std::vector<Parameters_> & params)
  : auxs_(0)
{
  ASSERT(odb.size() == params.size());
  for (std::size_t jobs = 0; jobs < params.size(); ++jobs) {
    auxs_.push_back(
      std::unique_ptr<ObsAuxControl_>(new ObsAuxControl_(odb[jobs], params[jobs])));
  }
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxControls<OBS>::ObsAuxControls(const ObsSpaces_ & odb, const eckit::Configuration & conf)
  : ObsAuxControls(odb,
                   // Split conf into subconfigurations, extract the "obs bias" section from each
                   // of them, then validate and deserialize that section into a Parameters_ object
                   validateAndDeserialize<Parameters_>(util::vectoriseAndFilter(conf, "obs bias")))
{}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxControls<OBS>::ObsAuxControls(const ObsAuxControls & other, const bool copy)
  : auxs_(other.size())
{
  Log::trace() << "ObsAuxControls<OBS>::ObsAuxControls copy starting" << std::endl;
  for (std::size_t jobs = 0; jobs < other.size(); ++jobs) {
    auxs_[jobs].reset(new ObsAuxControl_(other[jobs], copy));
  }
  Log::trace() << "ObsAuxControls<OBS>::ObsAuxControls copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxControls<OBS>::~ObsAuxControls() {
  Log::trace() << "ObsAuxControls<OBS>::~ObsAuxControls starting" << std::endl;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) auxs_[jobs].reset();
  Log::trace() << "ObsAuxControls<OBS>::~ObsAuxControls done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxControls<OBS>::read(const eckit::Configuration & conf) {
  Log::trace() << "ObsAuxControls<OBS>::read starting" << std::endl;
  std::vector<eckit::LocalConfiguration> obsconfs =
      conf.getSubConfigurations("observations.observers");
  ASSERT(obsconfs.size() == auxs_.size());
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    eckit::LocalConfiguration obsauxconf = obsconfs[jobs].getSubConfiguration("obs bias");
    Parameters_ params;
    params.validateAndDeserialize(obsauxconf);
    auxs_[jobs]->read(params);
  }
  Log::trace() << "ObsAuxControls<OBS>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxControls<OBS>::write(const eckit::Configuration & conf) const {
  Log::trace() << "ObsAuxControls<OBS>::write starting" << std::endl;
  std::vector<eckit::LocalConfiguration> obsconfs =
      conf.getSubConfigurations("observations.observers");
  ASSERT(obsconfs.size() == auxs_.size());
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    eckit::LocalConfiguration obsauxconf = obsconfs[jobs].getSubConfiguration("obs bias");
    Parameters_ params;
    params.validateAndDeserialize(obsauxconf);
    auxs_[jobs]->write(params);
  }
  Log::trace() << "ObsAuxControls<OBS>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
double ObsAuxControls<OBS>::norm() const {
  Log::trace() << "ObsAuxControls<OBS>::norm starting" << std::endl;
  double zz =  static_cast<double>(0.0);
  std::size_t ii = 0;
  double norm;
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    norm = auxs_[jobs]->norm();
    if (norm > 0.0) {
      zz += norm*norm;
      ++ii;
    }
  }
  Log::trace() << "ObsAuxControls<OBS>::norm done" << std::endl;
  return std::sqrt(zz/ii);
}

// -----------------------------------------------------------------------------

template<typename OBS>
ObsAuxControls<OBS> & ObsAuxControls<OBS>::operator=(const ObsAuxControls & rhs) {
  Log::trace() << "ObsAuxControl<OBS>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  for (std::size_t jobs = 0; jobs < auxs_.size(); ++jobs) {
    *auxs_[jobs] = *rhs.auxs_[jobs];
  }
  Log::trace() << "ObsAuxControl<OBS>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsAuxControls<OBS>::print(std::ostream & os) const {
  for (const auto & aux : auxs_) {
    os << *aux << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSAUXCONTROLS_H_
