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
#include "oops/interface/ModelData.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/Printable.h"

// -----------------------------------------------------------------------------

namespace detail {

oops::Variables modelVariables(const std::string modelName,
                               const oops::Variables & defaultModelVars,
                               const eckit::Configuration & config)
{
  oops::Variables returnVars;
  std::string includeVarsKey(modelName + " include variables");
  if (config.has(includeVarsKey)) {
    returnVars = oops::Variables(config.getStringVector(includeVarsKey));
  } else {
    std::string excludeVarsKey(modelName + " exclude variables");
    returnVars = defaultModelVars;
    if (config.has(excludeVarsKey)) {
        returnVars -= oops::Variables(config.getStringVector(excludeVarsKey));
    }
  }
  return returnVars;
}

}  // namespace detail

namespace oops {

// -----------------------------------------------------------------------------

/// Implementation of Geometry interface for a coupled model.
template <typename MODEL1, typename MODEL2>
class GeometryCoupled : public util::Printable {
 public:
  static std::string name1() {return MODEL1::name();}
  static std::string name2() {return MODEL2::name();}

  GeometryCoupled(const eckit::Configuration &, const eckit::mpi::Comm &);

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

  std::vector<size_t> variableSizes(const Variables & vars) const;

  const std::vector<Variables> & variables() const {return vars_;}

  bool levelsAreTopDown() const {return true;}
  const atlas::FunctionSpace & functionSpace() const {return nospace_;}
  const atlas::FieldSet & fields() const {return nofields_;}

 private:
  void print(std::ostream & os) const override;

  std::shared_ptr<Geometry<MODEL1>> geom1_;
  std::shared_ptr<Geometry<MODEL2>> geom2_;
  std::vector<Variables> vars_;  ///< variables that model1 and model2 should provide
  eckit::mpi::Comm * commPrints_;
  bool parallel_;
  int mymodel_;
  atlas::FunctionSpace nospace_;
  atlas::FieldSet nofields_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
GeometryCoupled<MODEL1, MODEL2>::GeometryCoupled(const eckit::Configuration & config,
                                                 const eckit::mpi::Comm & comm)
  : geom1_(), geom2_(), vars_(2),
    commPrints_(nullptr), parallel_(config.getBool("parallel", false)), mymodel_(-1)
{
  vars_[0] = ::detail::modelVariables(name1(),
                                      oops::ModelData<MODEL1>::defaultVariables(),
                                      config);
  vars_[1] = ::detail::modelVariables(name2(),
                                      oops::ModelData<MODEL2>::defaultVariables(),
                                      config);
  // check that the same variable isn't specified in both models' variables
  Variables commonvars = vars_[0];
  commonvars.intersection(vars_[1]);
  if (commonvars.size() > 0) {
    std::string errMsg = "Coupled model variable lists have overlap. "
                          "Use yaml to exclude these variables from one model:\n";
    for (auto variable : commonvars) {
        errMsg += variable.name() + "\n";
    }
    throw eckit::BadParameter(errMsg, Here());
  }
  if (parallel_) {
    const int mytask = comm.rank();
    const int ntasks = comm.size();
    const int tasks_per_model = ntasks / 2;
    mymodel_ = mytask / tasks_per_model + 1;

    // This creates the communicators for each model, named comm_model_{model name}
    // The first half of the MPI tasks will go to MODEL1, and the second half to MODEL2
    std::string commNameStr;
    if (mymodel_ == 1) commNameStr = "comm_model_" + name1();
    if (mymodel_ == 2) commNameStr = "comm_model_" + name2();
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commModel = comm.split(mymodel_, commName);

    if (mymodel_ == 1) {
      const eckit::LocalConfiguration conf1(config, name1());
      geom1_ = std::make_shared<Geometry<MODEL1>>(conf1, commModel);
    }
    if (mymodel_ == 2) {
      const eckit::LocalConfiguration conf2(config, name2());
      geom2_ = std::make_shared<Geometry<MODEL2>>(conf2, commModel);
    }

// This is creating Nprocs/2 new communicators, each of which pairs two processes:
// the N'th process among those handling model1 with the N'th process among
// those handling model2. This is used for handling prints.
    const int myrank = commModel.rank();

    std::string commPrintStr = "comm_ranks_" + std::to_string(myrank);
    char const *commPrintsName = commPrintStr.c_str();
    commPrints_ = &comm.split(myrank, commPrintsName);
  } else {
    const eckit::LocalConfiguration conf1(config, name1());
    geom1_ = std::make_shared<Geometry<MODEL1>>(conf1, comm);
    const eckit::LocalConfiguration conf2(config, name2());
    geom2_ = std::make_shared<Geometry<MODEL2>>(conf2, comm);
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
