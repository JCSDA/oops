/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_GEOMETRYFIXTURE_H_
#define TEST_INTERFACE_GEOMETRYFIXTURE_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

/// \brief Parameters loaded from the input YAML file and used by GeometryFixture.
template <typename MODEL>
class GeometryTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryTestParameters, Parameters)

 public:
  typedef typename oops::Geometry<MODEL>::Parameters_ GeometryParameters_;

  /// \brief Group of parameters controlling the tested model's geometry.
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};
  /// \brief Don't treat the presence of other parameter groups as an error (this makes it
  /// possible to reuse a single YAML file in tests of implementations of multiple oops interfaces).
  oops::IgnoreOtherParameters ignoreOthers{this};
};

// -----------------------------------------------------------------------------

/// \brief Fixture used by tests of the Geometry and GeometryIterator interfaces.
template <typename MODEL> class GeometryFixture : private boost::noncopyable {
  typedef oops::Geometry<MODEL>         Geometry_;
  typedef GeometryTestParameters<MODEL> TestParameters_;

 public:
  typedef typename Geometry_::Parameters_ Parameters_;
  static const Parameters_ & getParameters() {return getInstance().parameters_.geometry.value();}

 private:
  static GeometryFixture<MODEL>& getInstance() {
    static GeometryFixture<MODEL> theGeometryFixture;
    return theGeometryFixture;
  }

  GeometryFixture() {
    parameters_.validateAndDeserialize(TestEnvironment::config());
  }

  ~GeometryFixture() {}

  TestParameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_GEOMETRYFIXTURE_H_
