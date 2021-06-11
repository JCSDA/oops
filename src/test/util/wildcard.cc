/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/testing/Test.h"
#include "oops/util/wildcard.h"

namespace {

CASE("util/wildcard") {
  EXPECT(util::matchesWildcardPattern("", ""));
  EXPECT(util::matchesWildcardPattern("", "*"));
  EXPECT_NOT(util::matchesWildcardPattern("", "?"));
  EXPECT(util::matchesWildcardPattern("1", "1"));
  EXPECT(util::matchesWildcardPattern("12345", "12345"));
  EXPECT(util::matchesWildcardPattern("12345", "1?3?5"));
  EXPECT(util::matchesWildcardPattern("12345", "?2?4?"));
  EXPECT_NOT(util::matchesWildcardPattern("12345", "?3?4?"));
  EXPECT(util::matchesWildcardPattern("abcdefghijkl", "*cd*efg*k*"));
  EXPECT(util::matchesWildcardPattern("abcdefghijkl", "ab*cde*kl"));
  EXPECT_NOT(util::matchesWildcardPattern("abcdefghijkl", "*cd*defg*k*"));
  EXPECT_NOT(util::matchesWildcardPattern("abcdefghijkl", "ab*cde*k"));
  EXPECT_NOT(util::matchesWildcardPattern("abcdefghijkl", "b*cde*kl"));
  EXPECT(util::matchesWildcardPattern("abcdefghijkl", "***cd***efg***k***"));
  EXPECT(util::matchesWildcardPattern("abcdefghijkl", "ab***cde***kl"));
  EXPECT_NOT(util::matchesWildcardPattern("abcdefghijkl", "***cd***defg***k***"));
  EXPECT_NOT(util::matchesWildcardPattern("abcdefghijkl", "ab***cde***k"));
  EXPECT_NOT(util::matchesWildcardPattern("abcdefghijkl", "b***cde***kl"));
  EXPECT(util::matchesWildcardPattern("1212123", "*12123"));
  EXPECT(util::matchesWildcardPattern("1212123", "***12123"));
  EXPECT_NOT(util::matchesWildcardPattern("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "*a*a*a*b"));
  EXPECT_NOT(util::matchesWildcardPattern("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaab", "*a*a*a*a"));
  EXPECT(util::matchesWildcardPattern("xx", "*?"));
  EXPECT(util::matchesWildcardPattern("xx", "?*"));
}

}  // anonymous namespace

int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
