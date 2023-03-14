/*
 * (C) Copyright 2023- UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <vector>

#include "eckit/testing/Test.h"

#include "oops/base/Variables.h"
#include "oops/coupled/UtilsCoupled.h"
#include "oops/util/Logger.h"

namespace {

// -----------------------------------------------------------------------------

  CASE("test_split_existent") {
    using oops::Variables;
    using oops::splitVariables;
    typedef std::vector<Variables> VariablesVec;

    Variables vars1({"a", "b"});
    Variables vars2({"c", "d"});
    VariablesVec available({vars1, vars2});

    VariablesVec result;
    VariablesVec expected;

    result = splitVariables(Variables(), available);
    expected = VariablesVec({Variables(), Variables()});
    EXPECT(result == expected);

    result = splitVariables(Variables({"a"}), available);
    expected = VariablesVec({Variables({"a"}), Variables()});
    EXPECT(result == expected);

    result = splitVariables(Variables({"b", "a"}), available);
    expected = VariablesVec({Variables({"a", "b"}), Variables()});
    EXPECT(result == expected);

    result = splitVariables(Variables({"d"}), available);
    expected = VariablesVec({Variables(), Variables({"d"})});
    EXPECT(result == expected);

    result = splitVariables(Variables({"d", "a"}), available);
    expected = VariablesVec({Variables({"a"}), Variables({"d"})});
    EXPECT(result == expected);

    result = splitVariables(Variables({"d", "a", "b", "c"}), available);
    expected = VariablesVec({Variables({"a", "b"}), Variables({"c", "d"})});
    EXPECT(result == expected);
  }

// -----------------------------------------------------------------------------

  CASE("test_split_nonexistent") {
    using oops::Variables;
    using oops::splitVariables;
    typedef std::vector<Variables> VariablesVec;

    Variables vars1({"a", "b"});
    Variables vars2({"c", "d"});
    VariablesVec available({vars1, vars2});

    EXPECT_THROWS(splitVariables(Variables({"e"}), available));
    EXPECT_THROWS(splitVariables(Variables({"a", "d", "e"}), available));
    EXPECT_THROWS(splitVariables(Variables({"a", "b", "c", "d", "e"}), available));
  }

}  // anonymous namespace

// -----------------------------------------------------------------------------

int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
