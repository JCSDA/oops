/*
 * (C) Crown Copyright 2023, the Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/NormBase.h"
#include "oops/qg/QgTraits.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief Tests: Norm 1
/// Tests input of norm type.
///   i) Reads in norm type from configuration yaml file.
///  ii) Attempts to initialise norm object from factory.
/// iii) Norm type set in yaml file is intentionally invalid.
///  iv) Check that exception is raised for invalid norm type.
template <typename MODEL>
void testNorm1() {
  typedef oops::NormBase<MODEL> Norm_;

  const eckit::Configuration & config = TestEnvironment::config();

  std::unique_ptr<Norm_> norm_;

  const eckit::LocalConfiguration normConfig =
          config.getSubConfiguration("norm test 1").getSubConfiguration("norm");

  EXPECT_THROWS_MSG(norm_.reset(oops::NormFactory<MODEL>::create(normConfig)),
                    "Different does not exist in Norm factory.");
}

/// \brief Tests: Norm 2
/// Tests calculation of dot product using the norm multiplication methods.
///   i) Reads in geometry and norm type information from configuration yaml file.
///      Note: the value of alpha has been set to 2.0 in the yaml file.
///  ii) Creates two increments dx1 and dx2 from the geometry and sets all values to 1.0.
/// iii) Calculates the dot product of dx1.
///  iv) Performs multiplyMatrix and then multiplyMatrixInverse on dx1 and calculates dot product.
///   v) Checks that dot product of the new dx1 matches original dot product of dx1.
///  vi) Performs multiplyMatrixInverse and then multiplyMatrix on dx2 and calculates dot product.
/// vii) Checks that dot product of the new dx2 matches the new dot product of dx1 from (iv).
template <typename MODEL>
void testNorm2() {
    typedef oops::Increment<MODEL>   Increment_;
    typedef oops::Geometry<MODEL>    Geometry_;
    typedef oops::NormBase<MODEL>    Norm_;

    const eckit::Configuration & config = TestEnvironment::config();

    const Geometry_ geometry(config.getSubConfiguration("geometry"), oops::mpi::world());

    const oops::Variables incvars(config, "increment variables");
    const util::DateTime inctime(2020, 1, 1, 0, 0, 0);
    Increment_ dx1(geometry, incvars, inctime);
    dx1.random();
    Increment_ dx2(dx1);

    double dotProduct1 = dot_product(dx1, dx1);
    oops::Log::info() << "dot product 1 = " << dotProduct1 << std::endl;

    std::unique_ptr<Norm_> norm_;

    const eckit::LocalConfiguration normConfig =
            config.getSubConfiguration("norm test 2").getSubConfiguration("norm");

    ASSERT(normConfig.has("norm type"));
    norm_.reset(oops::NormFactory<MODEL>::create(normConfig));

    norm_->multiplyMatrix(dx1);
    norm_->multiplyMatrixInverse(dx1);

    double dotProduct2 = dot_product(dx1, dx1);
    oops::Log::info() << "dot product 2 = " << dotProduct2 << std::endl;

    const float tol = 1.0e-6;

    EXPECT(oops::is_close(static_cast<float>(dotProduct1), static_cast<float>(dotProduct2), tol));

    norm_->multiplyMatrixInverse(dx2);
    norm_->multiplyMatrix(dx2);

    double dotProduct3 = dot_product(dx2, dx2);
    oops::Log::info() << "dot product 3 = " << dotProduct3 << std::endl;

    EXPECT(oops::is_close(static_cast<float>(dotProduct2), static_cast<float>(dotProduct3), tol));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class TestNorm : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~TestNorm() = default;

 private:
  std::string testid() const override {return "test::TestNorm<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("assimilation/testNorm1")
      { testNorm1<MODEL>(); });

    ts.emplace_back(CASE("assimilation/testNorm2")
      { testNorm2<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
