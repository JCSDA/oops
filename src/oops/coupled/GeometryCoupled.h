/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Geometry.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// Parameters for Geometry describing a coupled model geometry
template<class... MODELs>
class GeometryCoupledParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryCoupledParameters, Parameters)

  /// Type of tuples stored in the GeometryCoupledParameters
  using RequiredParametersTupleT =
        std::tuple<RequiredParameter<typename Geometry<MODELs>::Parameters_>...>;
  /// Tuple that can be passed to the Parameter ctor
  using RequiredParameterInit = std::tuple<const char *, Parameters *>;

 public:
  /// Tuple of all Geometry Parameters.
  RequiredParametersTupleT geometries{RequiredParameterInit(MODELs::name().c_str(), this) ... };
  // Parameter to run the models sequentially or in parallel
  Parameter<bool> parallel{"parallel", false, this};
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

  /// Accessors to components of coupled geometry
  const Geometry<MODEL1> & geometry1() const {ASSERT(geom1_); return *geom1_;}
  const Geometry<MODEL2> & geometry2() const {ASSERT(geom2_); return *geom2_;}

  /// Accessors to model information
  const int & modelNumber() const {return mymodel_;}
  const bool & isParallel() const {return parallel_;}

  void latlon(std::vector<double> &, std::vector<double> &, const bool) const {}

 private:
  void print(std::ostream & os) const override;

  std::shared_ptr<Geometry<MODEL1>> geom1_;
  std::shared_ptr<Geometry<MODEL2>> geom2_;
  eckit::mpi::Comm * commPrints_;
  bool parallel_;
  int mymodel_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
GeometryCoupled<MODEL1, MODEL2>::GeometryCoupled(const Parameters_ & params,
                                                 const eckit::mpi::Comm & comm)
  : geom1_(), geom2_(), commPrints_(nullptr), parallel_(params.parallel.value()), mymodel_(-1)
{
  if (params.parallel) {
    Log::debug() << "Parallel coupled geometries" << std::endl;
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
      geom1_ = std::make_shared<Geometry<MODEL1>>(std::get<0>(params.geometries), commModel);
    }
    if (mymodel_ == 2) {
      geom2_ = std::make_shared<Geometry<MODEL2>>(std::get<1>(params.geometries), commModel);
    }

// This is creating Nprocs/2 new communicators, each of which pairs two processes:
// the N'th process among those handling model1 with the N'th process among
// those handling model2. This is used for handling prints.
    const int myrank = commModel.rank();

    std::string commPrintStr = "comm_ranks_" + std::to_string(myrank);
    char const *commPrintsName = commPrintStr.c_str();
    commPrints_ = &comm.split(myrank, commPrintsName);
  } else {
    Log::debug() << "Sequential coupled geometries" << std::endl;
    geom1_ = std::make_shared<Geometry<MODEL1>>(std::get<0>(params.geometries), comm);
    geom2_ = std::make_shared<Geometry<MODEL2>>(std::get<1>(params.geometries), comm);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void GeometryCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "GeometryCoupled::print starting" << std::endl;

  if (parallel_) {
    // Each model's geometry string is constructed on rank 0 for that model's communicator.
    // Here we gather the strings from each model onto the global rank 0 proc.
    std::stringstream ss;
    if (geom1_) {
      ss << std::endl << "GeometryCoupled: " << MODEL1::name() << std::endl;
      ss << std::setprecision(os.precision()) << *geom1_ << std::endl;
    }
    if (geom2_) {
      ss << std::endl << "GeometryCoupled: " << MODEL2::name() << std::endl;
      ss << std::setprecision(os.precision()) << *geom2_ << std::endl;
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
