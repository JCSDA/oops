/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "eckit/testing/Test.h"

#include "oops/util/DateTime.h"
#include "oops/util/PartialDateTime.h"


namespace {

  int notset = -1;

  CASE("test_toString_basic") {
    util::PartialDateTime pdt(2011, 9, 16, 13, 55, 20);
    std::string res = pdt.toString();
    std::string tar = "2011-09-16T13:55:20Z";
    EXPECT_EQUAL(res, tar);
  }

  CASE("test_toString_with_missing") {
    util::PartialDateTime pdt(2011, 9, 16, 13, -1, 20);
    std::string res = pdt.toString();
    std::string tar = "2011-09-16T13:**:20Z";
    EXPECT_EQUAL(res, tar);
  }

  CASE("test_construction_int_arg") {
    util::PartialDateTime pdt(2011, 9, 16, 13, 55, 20);
    EXPECT(pdt.year() == 2011);
    EXPECT(pdt.month() == 9);
    EXPECT(pdt.day() == 16);
    EXPECT(pdt.hour() == 13);
    EXPECT(pdt.minute() == 55);
    EXPECT(pdt.second() == 20);
  }


  CASE("test_construction_default") {
    util::PartialDateTime pdt{};
    EXPECT(pdt.year() == notset);
    EXPECT(pdt.month() == notset);
    EXPECT(pdt.day() == notset);
    EXPECT(pdt.hour() == notset);
    EXPECT(pdt.minute() == notset);
    EXPECT(pdt.second() == notset);
  }


  CASE("test_construction_string_arg") {
    util::PartialDateTime pdt("2011-09-16T13:55:20Z");
    EXPECT(pdt.year() == 2011);
    EXPECT(pdt.month() == 9);
    EXPECT(pdt.day() == 16);
    EXPECT(pdt.hour() == 13);
    EXPECT(pdt.minute() == 55);
    EXPECT(pdt.second() == 20);
  }


  CASE("test_construction_string_arg_unset") {
    // Check we properly handle the case where the string contains unset components
    // Note that any bit of a component results in it being considered unset.
    util::PartialDateTime pdt("2011-09-**T13:55:00Z");
    EXPECT(pdt.year() == 2011);
    EXPECT(pdt.month() == 9);
    EXPECT(pdt.day() == notset);
    EXPECT(pdt.hour() == 13);
    EXPECT(pdt.minute() == 55);
    EXPECT(pdt.second() == 0);

    pdt = util::PartialDateTime("2011-09-**T13:55:00Z");
    EXPECT(pdt.year() == 2011);
    EXPECT(pdt.month() == 9);
    EXPECT(pdt.day() == notset);
    EXPECT(pdt.hour() == 13);
    EXPECT(pdt.minute() == 55);
    EXPECT(pdt.second() == 0);
  }


  CASE("test_construction_string_bad_format") {
    // Check that we sanity check badly formatted strings

    // Defining a partially missing component
    EXPECT_THROWS(util::PartialDateTime("2011-09-*0T13:55:00Z"));

    // Wrong length
    EXPECT_THROWS(util::PartialDateTime("2011-09-01T13:55:00Z2"));
  }


  CASE("test_gt_operator") {
    // We also perform an equivelent check with DateTime operator working with a
    // PartialDateTime.
    util::DateTime dt1(2011, 9, 16, 13, 55, 20);

    util::PartialDateTime pdt1(notset, notset, notset, 12);
    util::PartialDateTime pdt2(notset, notset, notset, 13);
    util::PartialDateTime pdt3(notset, notset, notset, 14);
    util::PartialDateTime pdt4(2011, notset, notset, 14);
    util::PartialDateTime pdt5(2010, notset, notset, 14);

    EXPECT(!(pdt1 > dt1));
    EXPECT(!(dt1 < pdt1));

    EXPECT(!(pdt2 > dt1));
    EXPECT(!(dt1 < pdt2));

    EXPECT(pdt3 > dt1);
    EXPECT(dt1 < pdt3);

    EXPECT(pdt4 > dt1);
    EXPECT(dt1 < pdt4);

    EXPECT(!(pdt5 > dt1));
    EXPECT(!(dt1 < pdt5));
  }


  CASE("test_ge_operator") {
    // We also perform an equivelent check with DateTime operator working with a
    // PartialDateTime.
    util::DateTime dt1(2011, 9, 16, 13, 55, 20);

    util::PartialDateTime pdt1(notset, notset, notset, 12);
    util::PartialDateTime pdt2(notset, notset, notset, 13);
    util::PartialDateTime pdt3(notset, notset, notset, 14);
    util::PartialDateTime pdt4(2011, notset, notset, 14);
    util::PartialDateTime pdt5(2010, notset, notset, 14);

    EXPECT(!(pdt1 >= dt1));
    EXPECT(!(dt1 <= pdt1));

    EXPECT(pdt2 >= dt1);
    EXPECT(dt1 <= pdt2);

    EXPECT(pdt3 >= dt1);
    EXPECT(dt1 <= pdt3);

    EXPECT(pdt4 >= dt1);
    EXPECT(dt1 <= pdt4);

    EXPECT(!(pdt5 >= dt1));
    EXPECT(!(dt1 <= pdt5));
  }


  CASE("test_lt_operator") {
    // We also perform an equivelent check with DateTime operator working with a
    // PartialDateTime.
    util::DateTime dt1(2011, 9, 16, 13, 55, 20);

    util::PartialDateTime pdt1(notset, notset, notset, 12);
    util::PartialDateTime pdt2(notset, notset, notset, 13);
    util::PartialDateTime pdt3(notset, notset, notset, 14);
    util::PartialDateTime pdt4(2011, notset, notset, 14);
    util::PartialDateTime pdt5(2010, notset, notset, 14);

    EXPECT((pdt1 < dt1));
    EXPECT((dt1 > pdt1));

    EXPECT(!(pdt2 < dt1));
    EXPECT(!(dt1 > pdt2));

    EXPECT(!(pdt3 < dt1));
    EXPECT(!(dt1 > pdt3));

    EXPECT(!(pdt4 < dt1));
    EXPECT(!(dt1 > pdt4));

    EXPECT(pdt5 < dt1);
    EXPECT(dt1 > pdt5);
  }


  CASE("test_le_operator") {
    // We also perform an equivelent check with DateTime operator working with a
    // PartialDateTime.
    util::DateTime dt1(2011, 9, 16, 13, 55, 20);

    util::PartialDateTime pdt1(notset, notset, notset, 12);
    util::PartialDateTime pdt2(notset, notset, notset, 13);
    util::PartialDateTime pdt3(notset, notset, notset, 14);
    util::PartialDateTime pdt4(2011, notset, notset, 14);
    util::PartialDateTime pdt5(2010, notset, notset, 14);

    EXPECT((pdt1 <= dt1));
    EXPECT((dt1 >= pdt1));

    EXPECT(pdt2 <= dt1);
    EXPECT(dt1 >= pdt2);

    EXPECT(!(pdt3 <= dt1));
    EXPECT(!(dt1 >= pdt3));

    EXPECT(!(pdt4 <= dt1));
    EXPECT(!(dt1 >= pdt4));

    EXPECT(pdt5 <= dt1);
    EXPECT(dt1 >= pdt5);
  }


  CASE("test_eq_operator") {
    // We also perform an equivelent check with DateTime operator working with a
    // PartialDateTime.
    util::DateTime dt1(2011, 9, 16, 13, 55, 20);

    util::PartialDateTime pdt1(notset, notset, notset, 13);
    util::PartialDateTime pdt2(notset, notset, notset, 14);

    EXPECT(pdt1 == dt1);
    EXPECT(dt1 == pdt1);

    EXPECT(!(pdt2 == dt1));
    EXPECT(!(dt1 == pdt2));
  }


  CASE("test_ne_operator") {
    // We also perform an equivelent check with DateTime operator working with a
    // PartialDateTime.
    util::DateTime dt1(2011, 9, 16, 13, 55, 20);

    util::PartialDateTime pdt1(notset, notset, notset, 13);
    util::PartialDateTime pdt2(notset, notset, notset, 14);

    EXPECT(!(pdt1 != dt1));
    EXPECT(!(dt1 != pdt1));

    EXPECT((pdt2 != dt1));
    EXPECT((dt1 != pdt2));
  }


}  // namespace


int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
