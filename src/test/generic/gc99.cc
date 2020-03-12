/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cfloat>
#include <cmath>

#include "oops/generic/gc99.h"
#include "oops/util/Logger.h"

#include "eckit/testing/Test.h"

namespace {


// -----------------------------------------------------------------------------

  CASE("test_gc99") {
    double t;
    t = oops::gc99(0.0);
    oops::Log::info() << "gc99(0.0)=" << t << std::endl;
    EXPECT(std::abs(t-1) < 2*DBL_EPSILON);

    t = oops::gc99(0.5);
    oops::Log::info() << "gc99(0.5)=" << t << std::endl;
    EXPECT(std::abs(t-0.208333) < 1e-5);


    t = oops::gc99(1.0);
    oops::Log::info() << "gc99(1.0)=" << t << std::endl;
    EXPECT(std::abs(t) < 2*DBL_EPSILON);

    t = oops::gc99(2.0);
    oops::Log::info() << "gc99(2.0)=" << t << std::endl;
    EXPECT(std::abs(t) < 2*DBL_EPSILON);
  }

// -----------------------------------------------------------------------------

}  // anonymous namespace

int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
