/*
 * (C) Copyright 2023 UCAR-
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/FieldSet4D.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief Tests FieldSet3D ctor (from State and Increment fieldsets), valid time
/// and variables.
template <typename MODEL> void testFieldSet3D() {
  typedef oops::State<MODEL>     State_;
  typedef oops::Increment<MODEL> Increment_;
  typedef oops::Geometry<MODEL>  Geometry_;

  const eckit::Configuration & config = TestEnvironment::config();
  const bool parallel = config.getBool("parallel 4D");

  // Only run the 3D test in non-parallel mode.
  if (!parallel) {
    const Geometry_ geometry(config.getSubConfiguration("geometry"), oops::mpi::world());
    const State_ xx1(geometry, config.getSubConfiguration("state"));
    oops::Log::info() << "State: " << xx1 << std::endl;
    oops::FieldSet3D xx2Deep(xx1.fieldSet());
    oops::Log::info() << "FieldSet3D deep-copy: " << xx2Deep << std::endl;
    oops::FieldSet3D xx2Shallow(xx1.validTime(), xx1.geometry().getComm());
    xx2Shallow.shallowCopy(xx1.fieldSet());
    oops::Log::info() << "FieldSet3D shallow-copy: " << xx2Shallow << std::endl;

    // Check that the valid times are the same:
    EXPECT(xx2Deep.validTime() == xx1.validTime());
    EXPECT(xx2Shallow.validTime() == xx1.validTime());

    // Check that the variables are the same:
    oops::Log::info() << "State variables: " << xx1.variables() << std::endl;
    oops::Log::info() << "FieldSet3D variables: " << xx2Deep.variables() << std::endl;
    // Currently only comparing the variable strings since FieldSet3D provides levels in addition.
    // TODO(Algo): change to comparing .variables() when all models support State/Increment
    // variables that provide levels information.
    EXPECT(xx2Deep.variables().variables() == xx1.variables().variables());
    EXPECT(xx2Shallow.variables().variables() == xx1.variables().variables());

    const oops::Variables incvars(config, "increment variables");
    const util::DateTime inctime(2020, 1, 1, 0, 0, 0);
    Increment_ dx1(geometry, incvars, inctime);
    dx1.ones();
    oops::Log::info() << "Increment: " << dx1 << std::endl;
    const oops::FieldSet3D dx2(dx1.fieldSet());
    oops::Log::info() << "FieldSet3D: " << dx2 << std::endl;

    // Check that the valid times are the same:
    EXPECT(dx2.validTime() == dx1.validTime());

    // Check that the variables are the same:
    oops::Log::info() << "Increment variables: " << dx1.variables() << std::endl;
    oops::Log::info() << "FieldSet3D variables: " << dx2.variables() << std::endl;
    // Currently only comparing the variable strings since FieldSet3D provides levels in addition.
    // TODO(Algo): change to comparing .variables() when all models support State/Increment
    // variables that provide levels information.
    EXPECT(dx2.variables().variables() == dx1.variables().variables());

    // Test serialize-deserialize and initFieldSet3D
    std::vector<double> vect;
    dx2.serialize(vect);
    EXPECT(vect.size() > 0);
    EXPECT_EQUAL(vect.size(), dx2.serialSize());
    oops::FieldSet3D dx3 = oops::initFieldSet3D(dx2);
    dx3.zero();
    double dx2dp = dx2.dot_product_with(dx2, dx2.variables());
    double dx3dp = dx3.dot_product_with(dx3, dx3.variables());
    EXPECT_NOT_EQUAL(dx2dp, 0.0);
    EXPECT_EQUAL(dx3dp, 0.0);
    // Deserialize into dx3 and check that dx3 is now equal to dx2
    size_t index = 0;
    dx3.deserialize(vect, index);
    EXPECT_EQUAL(index, dx2.serialSize());
    EXPECT_EQUAL(index, dx3.serialSize());
    dx3dp = dx3.dot_product_with(dx3, dx3.variables());
    EXPECT_NOT_EQUAL(dx3dp, 0.0);
    EXPECT_EQUAL(dx3dp, dx3dp);
    // Check that an empty fieldset can't be deserialized into
    index = 0;
    oops::FieldSet3D empty(dx3.validTime(), dx3.commGeom());
    EXPECT_THROWS_AS(empty.deserialize(vect, index), eckit::BadParameter);
    // Do another deserialization
    dx2.serialize(vect);
    EXPECT_EQUAL(vect.size(), dx2.serialSize() * 2);
  }
}

// -----------------------------------------------------------------------------
/// \brief Tests FieldSet4D ctor (from State4D and Increment4D), valid times and
/// variables.
template <typename MODEL> void testFieldSet4D() {
  typedef oops::State4D<MODEL>     State4D_;
  typedef oops::Increment4D<MODEL> Increment4D_;
  typedef oops::Geometry<MODEL>    Geometry_;

  const eckit::Configuration & config = TestEnvironment::config();
  const bool parallel = config.getBool("parallel 4D");

  // Define space and time communicators
  const eckit::mpi::Comm * commSpace = &oops::mpi::world();
  const eckit::mpi::Comm * commTime = &oops::mpi::myself();
  if (parallel) {
    size_t ntasks = oops::mpi::world().size();
    size_t nslots = config.getInt("number of time slots per task");
    ASSERT(ntasks % nslots == 0);
    size_t myrank = oops::mpi::world().rank();
    size_t ntaskpslot = ntasks / nslots;
    size_t myslot = myrank / ntaskpslot;

    // Create a communicator for same sub-window, to be used for communications in space
    std::string sgeom = "comm_geom_" + std::to_string(myslot);
    char const *geomName = sgeom.c_str();
    commSpace = &oops::mpi::world().split(myslot, geomName);
    ASSERT(commSpace->size() == ntaskpslot);

    // Create a communicator for same local area, to be used for communications in time
    size_t myarea = commSpace->rank();
    std::string stime = "comm_time_" + std::to_string(myarea);
    char const *timeName = stime.c_str();
    commTime = &oops::mpi::world().split(myarea, timeName);
    ASSERT(commTime->size() == nslots);
  }

  const Geometry_ geometry(config.getSubConfiguration("geometry"), *commSpace, *commTime);
  const State4D_ xx1(geometry, config.getSubConfiguration("state4d"), *commTime);
  oops::Log::info() << "State4D: " << xx1 << std::endl;
  oops::FieldSet4D xx2(xx1);
  oops::Log::info() << "FieldSet4D: " << xx2 << std::endl;

  // Check that the valid times are the same:
  EXPECT(xx2.validTimes() == xx1.validTimes());

  // Check that the variables are the same:
  oops::Log::info() << "State4D variables: " << xx1.variables() << std::endl;
  oops::Log::info() << "FieldSet4D variables: " << xx2.variables() << std::endl;
  // Currently only comparing the variable strings since FieldSet3D provides levels in addition.
  // TODO(Algo): change to comparing .variables() when all models support State/Increment
  // variables that provide levels information.
  EXPECT(xx2.variables().variables() == xx1.variables().variables());

  // For Increment4D test, use zero increments with variables specified in yaml,
  // and the same times as in State4D
  const oops::Variables incvars(config, "increment variables");
  Increment4D_ dx1(geometry, incvars, xx1.times(), *commTime);
  dx1.ones();
  oops::Log::info() << "Increment4D: " << dx1 << std::endl;
  oops::FieldSet4D dx2(dx1);
  oops::Log::info() << "FieldSet4D: " << dx2 << std::endl;

  // Check that the valid times are the same:
  EXPECT(dx2.validTimes() == dx1.validTimes());

  // Check that the variables are the same:
  oops::Log::info() << "Increment4D variables: " << dx1.variables() << std::endl;
  oops::Log::info() << "FieldSet4D variables: " << dx2.variables() << std::endl;
  // Currently only comparing the variable strings since FieldSet3D provides levels in addition.
  // TODO(Algo): change to comparing .variables() when all models support State/Increment
  // variables that provide levels information.
  EXPECT(dx2.variables().variables() == dx1.variables().variables());
    // Test dot-product call
  const double norm1 = std::sqrt(dx1.dot_product_with(dx1));
  EXPECT(oops::is_close(norm1, dx2.norm(), 1.e-14));
  EXPECT_NOT_EQUAL(dx2.norm(), 0.0);
  // Test shallow and deep copies, +=, *= and zero method.
  // dx2_sc is a shallow copy of dx2; dx2_dc is a deep copy of dx2.
  oops::FieldSet4D dx2_sc = copyFieldSet4D(dx2, true);
  oops::FieldSet4D dx2_dc(dx2);
  EXPECT(oops::is_close(norm1, dx2_sc.norm(), 1.e-14));
  EXPECT(oops::is_close(norm1, dx2_dc.norm(), 1.e-14));
  // Set deep copy to zero; dx2 and dx2_sc should stay unchanged.
  dx2_dc.zero();
  EXPECT(oops::is_close(norm1, dx2_sc.norm(), 1.e-13));
  EXPECT_EQUAL(dx2_dc.norm(), 0.0);
  // Add zero fields to shallow copy and multiply by 1.0; should stay unchanged.
  dx2_sc += dx2_dc;
  dx2_sc *= 1.0;
  EXPECT(oops::is_close(norm1, dx2_sc.norm(), 1.e-13));
  // Multiply by zeroes shallow copy: both shallow copy and dx2 should now
  // be zero.
  dx2_sc *= dx2_dc;
  EXPECT_EQUAL(dx2.norm(), 0.0);
  EXPECT_EQUAL(dx2_sc.norm(), 0.0);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class FieldSet4D : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~FieldSet4D() = default;

 private:
  std::string testid() const override {return "test::FieldSet4D<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/FieldSet4D/testFieldSet3D")
      { testFieldSet3D<MODEL>(); });

    ts.emplace_back(CASE("interface/FieldSet4D/testFieldSet4D")
      { testFieldSet4D<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
