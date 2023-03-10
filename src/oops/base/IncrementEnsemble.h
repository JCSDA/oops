/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENTENSEMBLE_H_
#define OOPS_BASE_INCREMENTENSEMBLE_H_

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template <typename MODEL> class GeometryIterator;

/// Parameters for ensemble members from template.
template <typename MODEL>
class IncrementMemberTemplateParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementMemberTemplateParameters, Parameters)

  typedef typename Increment<MODEL>::ReadParameters_ Parameters_;
 public:
  RequiredParameter<Parameters_> increment{"template", "template to define a generic member", this};
  RequiredParameter<std::string> pattern{"pattern", "pattern to be replaced for members", this};
  RequiredParameter<size_t> nmembers{"nmembers", "number of members", this};
  Parameter<size_t> start{"start", "starting member index", 1, this};
  Parameter<std::vector<size_t>> except{"except", "excluded members indices", {}, this};
  Parameter<size_t> zpad{"zero padding", "zero padding", 0, this};
};

// -----------------------------------------------------------------------------

/// Parameters for the ensemble of increments.
template <typename MODEL>
class IncrementEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementEnsembleParameters, Parameters)

  typedef typename Increment<MODEL>::ReadParameters_ Parameters_;
  typedef IncrementMemberTemplateParameters<MODEL> IncrementMemberTemplateParameters_;
 public:
  RequiredParameter<util::DateTime> date{"date", "increment date", this};
  OptionalParameter<std::vector<Parameters_>> increments{"members",
                   "members of the increment ensemble", this};
  OptionalParameter<IncrementMemberTemplateParameters_> increments_template{"members from template",
                   "template to define members of the increment ensemble", this};

  /// Overridden to detect missing conditionally required parameters
  using Parameters::deserialize;
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  /// Get ensemble size
  size_t size() const;

  /// Get Increment parameters for a given ensemble index
  Parameters_ getIncrementParameters(const size_t &) const;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementEnsembleParameters<MODEL>::deserialize(util::CompositePath &path,
                                                     const eckit::Configuration &config)
{
  Parameters::deserialize(path, config);

  if (increments.value() == boost::none && increments_template.value() == boost::none) {
    throw eckit::UserError(
        path.path() +
        ": both members and members from template are missing",
        Here());
  }
  if (increments.value() != boost::none && increments_template.value() != boost::none) {
    throw eckit::UserError(
        path.path() +
        ": both members and members from template are present",
        Here());
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
size_t IncrementEnsembleParameters<MODEL>::size() const
{
  if (increments.value() != boost::none) {
    return increments.value()->size();
  } else {
    return increments_template.value()->nmembers.value();
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
typename Increment<MODEL>::ReadParameters_
  IncrementEnsembleParameters<MODEL>::getIncrementParameters(const size_t & ie) const
{
  // Check ensemble size
  if (ie >= this->size()) {
    ABORT("IncrementEnsembleParameters: getIncrementParameters member index is too large");
  }

  if (increments.value() != boost::none) {
    // Explicit members
    return (*increments.value())[ie];
  } else {
    // Members template

    // Template configuration
    eckit::LocalConfiguration incrementConf;
    increments_template.value()->increment.value().serialize(incrementConf);

    // Get correct index
    size_t count = increments_template.value()->start;
    for (size_t jj = 0; jj <= ie; ++jj) {
      // Check for excluded members
      while (std::count(increments_template.value()->except.value().begin(),
             increments_template.value()->except.value().end(), count)) {
        count += 1;
      }

      // Update counter
      if (jj < ie) count += 1;
    }

    // Copy and update template configuration with pattern
    eckit::LocalConfiguration memberConf(incrementConf);

    // Replace pattern recursively in the configuration
    util::seekAndReplace(memberConf, increments_template.value()->pattern,
      count, increments_template.value()->zpad);

    // Get member parameters
    Parameters_ params;
    params.validateAndDeserialize(memberConf);
    return params;
  }
}

// -----------------------------------------------------------------------------
/// Parameters for the ensemble of increments generated from ensemble of states
/// with specified inflation and linear variable changes.
template <typename MODEL>
class IncrementEnsembleFromStatesParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementEnsembleFromStatesParameters, Parameters)

  typedef typename Increment<MODEL>::ReadParameters_ IncrementReadParameters_;
  typedef typename LinearVariableChange<MODEL>::Parameters_ LinearVarChangeParameters_;
 public:
  OptionalParameter<IncrementReadParameters_> inflationField{"inflation field",
                   "inflation field (as increment in model space)", this};
  Parameter<double> inflationValue{"inflation value", "inflation value (scalar)",
                    1.0, this};
  OptionalParameter<LinearVarChangeParameters_> linVarChange{"linear variable change",
                   "linear variable changes applied to the increments", this};
  StateEnsembleParameters<MODEL> states{this};
};


/// \brief Ensemble of increments
template<typename MODEL> class IncrementEnsemble {
  typedef Geometry<MODEL>                    Geometry_;
  typedef GeometryIterator<MODEL>            GeometryIterator_;
  typedef Increment<MODEL>                   Increment_;
  typedef LinearVariableChange<MODEL>        LinearVariableChange_;
  typedef State<MODEL>                       State_;
  typedef StateEnsemble<MODEL>               StateEnsemble_;
  typedef IncrementEnsembleFromStatesParameters<MODEL> IncrementEnsembleFromStatesParameters_;
  typedef IncrementEnsembleParameters<MODEL> IncrementEnsembleParameters_;
  typedef StateEnsembleParameters<MODEL>     StateEnsembleParameters_;
  typedef StateParametersND<MODEL>           StateParameters_;

 public:
  /// Constructor
  IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                    const util::DateTime &, const int rank);
  IncrementEnsemble(const IncrementEnsembleFromStatesParameters_ &, const State_ &, const State_ &,
                    const Geometry_ &, const Variables &);
  /// \brief construct ensemble of perturbations by reading them from disk
  IncrementEnsemble(const Geometry_ &, const Variables &, const IncrementEnsembleParameters_ &);
  /// \brief construct ensemble of perturbations by reading two state ensembles (one member at a
  //         time) and taking the  difference of each set of pairs
  IncrementEnsemble(const Geometry_ &, const Variables &, const StateEnsembleParameters_ &,
                    const StateEnsembleParameters_ &);

  void write(const eckit::Configuration &) const;

  /// Accessors
  size_t size() const {return ensemblePerturbs_.size();}
  Increment_ & operator[](const int ii) {return ensemblePerturbs_[ii];}
  const Increment_ & operator[](const int ii) const {return ensemblePerturbs_[ii];}

  void packEigen(Eigen::MatrixXd & X, const GeometryIterator_ & gi) const;
  void setEigen(const Eigen::MatrixXd & X, const GeometryIterator_ & gi);

 private:
  std::vector<Increment_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                                            const util::DateTime & tslot, const int rank)
  : ensemblePerturbs_()
{
  ensemblePerturbs_.reserve(rank);
  for (int m = 0; m < rank; ++m) {
    ensemblePerturbs_.emplace_back(resol, vars, tslot);
  }
  Log::trace() << "IncrementEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const IncrementEnsembleFromStatesParameters_ & params,
                                            const State_ & xb, const State_ & fg,
                                            const Geometry_ & resol, const Variables & vars)
  : ensemblePerturbs_()
{
  Log::trace() << "IncrementEnsemble:contructor start" << std::endl;
  // Check sizes and fill in timeslots
  util::DateTime tslot = xb.validTime();

  // Read inflation field
  std::unique_ptr<Increment_> inflationField;
  if (params.inflationField.value() != boost::none) {
    inflationField.reset(new Increment_(resol, vars, tslot));
    inflationField->read(*params.inflationField.value());
  }

  // Get inflation value
  const double inflationValue = params.inflationValue;

  // Setup change of variable
  std::unique_ptr<LinearVariableChange_> linvarchg;
  if (params.linVarChange.value() != boost::none) {
    const auto & linvar = *params.linVarChange.value();
    linvarchg.reset(new LinearVariableChange_(resol, linvar));
    linvarchg->changeVarTraj(fg, *linvar.outputVariables.value());
  }

  // Create a zero state (to be kept constant for building increments)
  State_ zerov = State_(resol, vars, tslot);
  zerov.zero();

  // Initialize the ensemble mean accumulator
  Accumulator<MODEL, State_, State_> ensmean(zerov);

  // Ensemble size from parameters and its inverse for the mean
  const size_t nens = params.states.size();
  const double rr = 1.0/static_cast<double>(nens);

  // Reserve memory for the departures ensemble
  ensemblePerturbs_.reserve(nens);
  for (size_t jj = 0; jj < nens; ++jj) {
    Increment_ dx(resol, vars, tslot);
    ensemblePerturbs_.emplace_back(std::move(dx));
  }

  // Read the state, compute the mean and store states at increment resolution.
  for (size_t jj = 0; jj < nens; ++jj) {
    // 1. Read state jj
    State_ xx(resol, params.states.getStateParameters(jj));
    // 2. Accumulate it to mean.
    ensmean.accumul(rr, xx);
    // 3. Subtract zerov so to have an increment at ir res. and push to ensPert
    ensemblePerturbs_[jj].diff(xx, zerov);
  }
  State_ bgmean = std::move(ensmean);

  // 4. Subtract zerov from bgmean so to have an increment at ir res. and store in bgmean_ir
  Increment_ bgmean_ir(resol, vars, tslot);
  bgmean_ir.diff(bgmean, zerov);

  // 5. Loop on the stored states at increment resolution
  for (size_t jj = 0; jj < nens; ++jj) {
    // 5.1 remove the mean
    ensemblePerturbs_[jj] -= bgmean_ir;
    // 5.2 Apply inflation
    if (params.inflationField.value() != boost::none) {
      ensemblePerturbs_[jj].schur_product_with(*inflationField);
    }
    ensemblePerturbs_[jj] *= inflationValue;
    // 5.3 Apply lin chvar
    if (params.linVarChange.value() != boost::none) {
      const auto & linvar = *params.linVarChange.value();
      oops::Variables varin = *linvar.inputVariables.value();
      linvarchg->changeVarInverseTL(ensemblePerturbs_[jj], varin);
    }
  }
  Log::trace() << "IncrementEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                                            const IncrementEnsembleParameters_ & params)
  : ensemblePerturbs_()
{
  // Datetime for ensemble
  util::DateTime tslot = params.date;

  // Reserve memory to hold ensemble
  const size_t nens = params.size();
  ensemblePerturbs_.reserve(nens);

  // Loop over all ensemble members
  for (size_t jj = 0; jj < nens; ++jj) {
    Increment_ dx(resol, vars, tslot);
    dx.read(params.getIncrementParameters(jj));
    ensemblePerturbs_.emplace_back(std::move(dx));
  }
  Log::trace() << "IncrementEnsemble:contructor (by reading increment ensemble) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                                            const StateEnsembleParameters_ & configBase,
                                            const StateEnsembleParameters_ & configPert)
  : ensemblePerturbs_()
{
  // Check base and perturbations ensembles size
  size_t baseNens = configBase.size();
  size_t pertNens = configPert.size();
  ASSERT(baseNens == pertNens);

  // Reserve memory to hold ensemble
  ensemblePerturbs_.reserve(baseNens);

  // Loop over all ensemble members
  for (size_t jj = 0; jj < baseNens; ++jj) {
    // Load base member
    State_ xBase(resol, configBase.getStateParameters(jj));

    // Load perturbation member
    State_ xPert(resol, configPert.getStateParameters(jj));

    // Difference
    Increment_ dx(resol, vars, xBase.validTime());
    dx.diff(xBase, xPert);
    ensemblePerturbs_.emplace_back(std::move(dx));
  }
  Log::trace() << "IncrementEnsemble:contructor (by diffing state ensembles) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble<MODEL>::write(const eckit::Configuration & config) const
{
  eckit::LocalConfiguration outConfig(config);
  for (size_t ii=0; ii < size(); ++ii) {
    outConfig.set("member", ii+1);
    ensemblePerturbs_[ii].write(outConfig);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble<MODEL>::packEigen(Eigen::MatrixXd & X,
                                         const GeometryIterator_ & gi) const
{
  size_t ngp = ensemblePerturbs_[0].getLocal(gi).getVals().size();
  size_t nens = ensemblePerturbs_.size();
  X.resize(ngp, nens);
  for (size_t jj=0; jj < nens; ++jj) {
    LocalIncrement gp = ensemblePerturbs_[jj].getLocal(gi);
    std::vector<double> tmp1 = gp.getVals();
    for (size_t iv=0; iv < ngp; ++iv) {
      X(iv, jj) = tmp1[iv];
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble<MODEL>::setEigen(const Eigen::MatrixXd & X,
                                        const GeometryIterator_ & gi)
{
  size_t ngp = ensemblePerturbs_[0].getLocal(gi).getVals().size();
  size_t nens = ensemblePerturbs_.size();

  LocalIncrement gptmp = ensemblePerturbs_[0].getLocal(gi);
  std::vector<double> tmp = gptmp.getVals();

  for (size_t jj=0; jj < nens; ++jj) {
    for (size_t iv=0; iv < ngp; ++iv) {
      tmp[iv] = X(iv, jj);
    }
    gptmp.setVals(tmp);
    ensemblePerturbs_[jj].setLocal(gptmp, gi);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INCREMENTENSEMBLE_H_
