/*
 * (C) Copyright 2021-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/coupled/UtilsCoupled.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Parameters for Geometry describing a coupled model geometry
template<typename MODEL1, typename MODEL2>
class GeometryCoupledParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryCoupledParameters, Parameters)

  typedef typename Geometry<MODEL1>::Parameters_ Parameters1_;
  typedef typename Geometry<MODEL2>::Parameters_ Parameters2_;

 public:
  RequiredParameter<Parameters1_> geometry1{MODEL1::name().c_str(), this};
  RequiredParameter<Parameters2_> geometry2{MODEL2::name().c_str(), this};
  Parameter<bool> parallel{"parallel", "run the models sequentially or in parallel",
                           false, this};
  Parameter<Variables> vars1{std::string(MODEL1::name() + " variables").c_str(),
          "variables that the first model should provide, have to be different "
          "from the variables that the second model provides", {}, this};
  Parameter<Variables> vars2{std::string(MODEL2::name() + " variables").c_str(),
          "variables that the second model should provide, have to be different "
          "from the variables that the first model provides", {}, this};
};

// -----------------------------------------------------------------------------

/// Implementation of Geometry interface for a coupled model.
template <typename MODEL1, typename MODEL2>
class GeometryCoupled : public util::Printable {
 public:
  typedef GeometryCoupledParameters<MODEL1, MODEL2> Parameters_;

  GeometryCoupled(const Parameters_ &, const eckit::mpi::Comm &);

  /// Accessor to the MPI communicator between models
  const eckit::mpi::Comm & getCommPairRanks() const {ASSERT(commPrints_); return *commPrints_;}

  /// WARNING: This implementation is wrong because there are in general two communicators.
  ///          It is provided for compile-time compatibility with oops interfaces, but will throw
  ///          an exception if called as a reminder that the implementation is incorrect.
  const eckit::mpi::Comm & getComm() const {
    throw eckit::Exception("Called GeometryCoupled.getComm(), but this is just a stub");
    return oops::mpi::world();
  }

  /// Accessors to components of coupled geometry
  const Geometry<MODEL1> & geometry1() const {ASSERT(geom1_); return *geom1_;}
  const Geometry<MODEL2> & geometry2() const {ASSERT(geom2_); return *geom2_;}

  /// Accessors to model information
  const int & modelNumber() const {return mymodel_;}
  const bool & isParallel() const {return parallel_;}

  void latlon(std::vector<double> &, std::vector<double> &, const bool) const {}

  std::vector<size_t> variableSizes(const Variables & vars) const;

  const std::vector<Variables> & variables() const {return vars_;}

  bool levelsAreTopDown() const {return true;}
  const atlas::FunctionSpace & functionSpace() const {return nospace_;}
  const atlas::FieldSet & extraFields() const {return nofields_;}

 private:
  void print(std::ostream & os) const override;

  std::shared_ptr<Geometry<MODEL1>> geom1_;
  std::shared_ptr<Geometry<MODEL2>> geom2_;
  const std::vector<Variables> vars_;  ///< variables that model1 and model2 should provide
  eckit::mpi::Comm * commPrints_;
  bool parallel_;
  int mymodel_;
  atlas::FunctionSpace nospace_;
  atlas::FieldSet nofields_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
GeometryCoupled<MODEL1, MODEL2>::GeometryCoupled(const Parameters_ & params,
                                                 const eckit::mpi::Comm & comm)
  : geom1_(), geom2_(), vars_({params.vars1.value(), params.vars2.value()}),
    commPrints_(nullptr), parallel_(params.parallel.value()), mymodel_(-1)
{
  // check that the same variable isn't specified in both models'
  // variables
  Variables commonvars = vars_[0];
  commonvars.intersection(vars_[1]);
  if (commonvars.size() > 0) {
    throw eckit::BadParameter("Variables for different components of coupled "
          "model can not overlap", Here());
  }
  if (params.parallel) {
    const int mytask = comm.rank();
    const int ntasks = comm.size();
    const int tasks_per_model = ntasks / 2;
    mymodel_ = mytask / tasks_per_model + 1;

    // This creates the communicators for each model, named comm_model_{model name}
    // The first half of the MPI tasks will go to MODEL1, and the second half to MODEL2
    std::string commNameStr;
    if (mymodel_ == 1) commNameStr = "comm_model_" + MODEL1::name();
    if (mymodel_ == 2) commNameStr = "comm_model_" + MODEL2::name();
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commModel = comm.split(mymodel_, commName);

    if (mymodel_ == 1) {
      geom1_ = std::make_shared<Geometry<MODEL1>>(params.geometry1, commModel);
    }
    if (mymodel_ == 2) {
      geom2_ = std::make_shared<Geometry<MODEL2>>(params.geometry2, commModel);
    }

// This is creating Nprocs/2 new communicators, each of which pairs two processes:
// the N'th process among those handling model1 with the N'th process among
// those handling model2. This is used for handling prints.
    const int myrank = commModel.rank();

    std::string commPrintStr = "comm_ranks_" + std::to_string(myrank);
    char const *commPrintsName = commPrintStr.c_str();
    commPrints_ = &comm.split(myrank, commPrintsName);
  } else {
    geom1_ = std::make_shared<Geometry<MODEL1>>(params.geometry1, comm);
    geom2_ = std::make_shared<Geometry<MODEL2>>(params.geometry2, comm);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
std::vector<size_t> GeometryCoupled<MODEL1, MODEL2>::variableSizes(const Variables & vars) const {
  // decide what variables are provided by what model
  std::vector<Variables> splitvars = splitVariables(vars, vars_);
  const std::vector<size_t> reqvars1sizes = geom1_->variableSizes(splitvars[0]);
  const std::vector<size_t> reqvars2sizes = geom2_->variableSizes(splitvars[1]);

  std::vector<size_t> varsizes(vars.size());
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    if (splitvars[0].has(vars[jvar])) {
      varsizes[jvar] = reqvars1sizes[splitvars[0].find(vars[jvar])];
    } else {
      varsizes[jvar] = reqvars2sizes[splitvars[1].find(vars[jvar])];
    }
  }
  return varsizes;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void GeometryCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "GeometryCoupled::print starting" << std::endl;

  if (parallel_) {
    // Each model's geometry string is constructed on rank 0 for that model's communicator.
    // Here we gather the strings from each model onto the global rank 0 proc.
    std::stringstream ss;
    ss.setf(os.flags());
    ss.precision(os.precision());
    if (geom1_) {
      ss << std::endl << "GeometryCoupled: " << MODEL1::name() << std::endl;
      ss << *geom1_ << std::endl;
    }
    if (geom2_) {
      ss << std::endl << "GeometryCoupled: " << MODEL2::name() << std::endl;
      ss << *geom2_ << std::endl;
    }
    util::gatherPrint(os, ss.str(), *commPrints_);
  } else {
    os << std::endl << "GeometryCoupled: " << MODEL1::name() << std::endl;
    os << *geom1_ << std::endl;
    os << std::endl << "GeometryCoupled: " << MODEL2::name() << std::endl;
    os << *geom2_;
  }

  Log::trace() << "GeometryCoupled::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
