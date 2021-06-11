/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cfloat>
#include <cmath>

#include "oops/generic/soar.h"
#include "oops/util/Logger.h"

#include "eckit/testing/Test.h"

namespace {


// -----------------------------------------------------------------------------

  CASE("test_soar") {
    double t;
    t = oops::soar(0.0);
    oops::Log::info() << "soar(0.0)=" << t << std::endl;
    EXPECT(std::abs(t-1) < 2*DBL_EPSILON);

    t = oops::soar(1.0);
    oops::Log::info() << "soar(1.0)=" << t << std::endl;
    EXPECT(std::abs(t-0.735758882342885) < 1e-6);

    t = oops::soar(2.0);
    oops::Log::info() << "soar(2.0)=" << t << std::endl;
    EXPECT(std::abs(t-0.406005849709838) < 1e-6);

    t = oops::soar(4.0);
    oops::Log::info() << "soar(4.0)=" << t << std::endl;
    EXPECT(std::abs(t-0.091578194443671) < 1e-6);

    t = oops::soar(8.0);
    oops::Log::info() << "soar(8.0)=" << t << std::endl;
    EXPECT(std::abs(t-0.003019163651123) < 1e-6);
  }

// -----------------------------------------------------------------------------

}  // anonymous namespace

int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
