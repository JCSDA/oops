/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <tuple>
#include <utility>

#include "eckit/mpi/Comm.h"

#include "oops/base/Geometry.h"
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
};

// -----------------------------------------------------------------------------

/// Implementation of Geometry interface for a coupled model.
template <typename MODEL1, typename MODEL2>
class GeometryCoupled : public util::Printable {
 public:
  typedef GeometryCoupledParameters<MODEL1, MODEL2> Parameters_;

  GeometryCoupled(const Parameters_ &, const eckit::mpi::Comm &);

  /// Accessor to the MPI communicator
  const eckit::mpi::Comm & getComm() const {return comm_;}

  /// Accessors to components of coupled geometry
  const Geometry<MODEL1> & geometry1() const {return *geom1_;}
  const Geometry<MODEL2> & geometry2() const {return *geom2_;}

 private:
  void print(std::ostream & os) const override;

  std::shared_ptr<Geometry<MODEL1>> geom1_;
  std::shared_ptr<Geometry<MODEL2>> geom2_;
  const eckit::mpi::Comm & comm_;
};

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2>
GeometryCoupled<MODEL1, MODEL2>::GeometryCoupled(const Parameters_ & params,
                                                 const eckit::mpi::Comm & comm)
  : geom1_(), geom2_(), comm_(comm)
{
  geom1_ = std::make_shared<Geometry<MODEL1>>(std::get<0>(params.geometries), comm);
  geom2_ = std::make_shared<Geometry<MODEL2>>(std::get<1>(params.geometries), comm);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void GeometryCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  os << "GeometryCoupled: " << MODEL1::name() << std::endl;
  os << *geom1_ << std::endl;
  os << "GeometryCoupled: " << MODEL2::name() << std::endl;
  os << *geom2_;
}

// -----------------------------------------------------------------------------

}  // namespace oops
