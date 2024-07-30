/*
 * (C) Crown copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_TIMEWINDOW_H_
#define TEST_UTIL_TIMEWINDOW_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/TimeWindow.h"

namespace test {

  void testConstructor() {
    // Set up time window using "begin", "length" and the default value of "bound to include".
    const util::DateTime windowStart(2023, 1, 1, 0, 0, 0);
    const util::Duration windowLength("P1D");
    eckit::LocalConfiguration configBeginLength;
    configBeginLength.set("begin", windowStart.toString());
    configBeginLength.set("length", windowLength.toString());
    const util::TimeWindow timeWindowBeginLength(configBeginLength);

    // Set up time window using "begin", "length" and "bound to include" set to "begin".
    configBeginLength.set("bound to include", "begin");
    const util::TimeWindow timeWindowBeginLengthLower(configBeginLength);

    // Set up time window using "begin", "length" and "bound to include" set to "end".
    configBeginLength.set("bound to include", "end");
    const util::TimeWindow timeWindowBeginLengthUpper(configBeginLength);

    // Set up time window using "begin", "length" and an invalid value of "bound to include".
    eckit::LocalConfiguration configBeginLengthInvalid;
    configBeginLengthInvalid.set("begin", windowStart.toString());
    configBeginLengthInvalid.set("length", windowLength.toString());
    configBeginLengthInvalid.set("bound to include", "invalid");
    EXPECT_THROWS(const util::TimeWindow timeWindowBeginLengthInvalid(configBeginLengthInvalid));

    // Set up time window using "begin", "end" and the default value of "bound to include".
    const util::DateTime windowEnd(2021, 2, 1, 0, 0, 0);
    eckit::LocalConfiguration configBeginEnd;
    configBeginEnd.set("begin", windowStart.toString());
    configBeginEnd.set("end", windowEnd.toString());
    const util::TimeWindow timeWindowBeginEnd(configBeginEnd);

    // Set up time window using "begin", "end" and "bound to include" set to "begin".
    configBeginEnd.set("bound to include", "begin");
    const util::TimeWindow timeWindowBeginEndLower(configBeginEnd);

    // Set up time window using "begin", "end" and "bound to include" set to "end".
    configBeginEnd.set("bound to include", "end");
    const util::TimeWindow timeWindowBeginEndUpper(configBeginEnd);

    // Set up time window using "begin", "length" and "end", throwing an exception.
    eckit::LocalConfiguration configBeginLengthEnd;
    configBeginLengthEnd.set("begin", windowStart.toString());
    configBeginLengthEnd.set("length", windowLength.toString());
    configBeginLengthEnd.set("end", windowEnd.toString());
    EXPECT_THROWS(const util::TimeWindow timeWindowBeginLengthEnd(configBeginLengthEnd));
  }

  void testAccessors() {
    // Set up time window.
    const util::DateTime windowStart(2023, 1, 1, 0, 0, 0);
    const util::Duration windowLength("P1D");
    eckit::LocalConfiguration config;
    config.set("begin", windowStart.toString());
    config.set("length", windowLength.toString());
    const util::TimeWindow timeWindow(config);

    // Test routines that access information about the time window.
    const util::DateTime expectedWindowStart(windowStart);
    const util::Duration expectedWindowLength(windowLength);
    const util::DateTime expectedWindowEnd(expectedWindowStart + expectedWindowLength);
    const util::DateTime expectedWindowMidpoint
      (expectedWindowStart + expectedWindowLength / 2);

    EXPECT(timeWindow.start() == expectedWindowStart);
    EXPECT(timeWindow.end() == expectedWindowEnd);
    EXPECT(timeWindow.length() == expectedWindowLength);
    EXPECT(timeWindow.midpoint() == expectedWindowMidpoint);

    // Configure the TimeWindow object to use an inclusive lower bound and repeat the tests.
    // Expect no difference in the results.
    config.set("bound to include", "begin");
    const util::TimeWindow timeWindowInclusiveLower(config);
    EXPECT(timeWindowInclusiveLower.start() == expectedWindowStart);
    EXPECT(timeWindowInclusiveLower.end() == expectedWindowEnd);
    EXPECT(timeWindowInclusiveLower.length() == expectedWindowLength);
    EXPECT(timeWindowInclusiveLower.midpoint() == expectedWindowMidpoint);
  }

  void testSubWindow() {
    // Set up time window.
    const util::DateTime windowStart(2023, 1, 1, 0, 0, 0);
    const util::Duration windowLength("P1D");
    eckit::LocalConfiguration config;
    config.set("begin", windowStart.toString());
    config.set("length", windowLength.toString());

    const util::TimeWindow timeWindow(config);

    // Sub-window midpoint.
    const util::DateTime midPoint(2023, 1, 1, 12, 0, 0);
    // Window around midpoint in which observations are accepted.
    const util::Duration halfWidth("PT3H");
    // Create sub-window spanning `midPoint` +/- `halfHidth`.
    const util::TimeWindow subWindow = timeWindow.createSubWindow(midPoint, halfWidth);

    // Check sub-window accessors.
    const util::DateTime expectedSubWindowStart(midPoint - halfWidth);
    const util::DateTime expectedSubWindowEnd(midPoint + halfWidth);
    const util::Duration expectedSubWindowLength(expectedSubWindowEnd - expectedSubWindowStart);
    const util::DateTime expectedSubWindowMidpoint
      (expectedSubWindowStart + expectedSubWindowLength / 2);

    EXPECT(subWindow.start() == expectedSubWindowStart);
    EXPECT(subWindow.end() == expectedSubWindowEnd);
    EXPECT(subWindow.length() == expectedSubWindowLength);
    EXPECT(subWindow.midpoint() == expectedSubWindowMidpoint);

    // Perform the same checks with a large half-width. In this case expect the
    // sub-window to have the same bounds as the original window.
    const util::Duration halfWidthLarge("P1D");
    const util::TimeWindow subWindowLarge = timeWindow.createSubWindow(midPoint, halfWidthLarge);

    EXPECT(subWindowLarge.start() == timeWindow.start());
    EXPECT(subWindowLarge.end() == timeWindow.end());
    EXPECT(subWindowLarge.length() == timeWindow.length());
    EXPECT(subWindowLarge.midpoint() == timeWindow.midpoint());

    // Create another sub-window spanning `subBegin`, `subBegin` + `subLength`.
    const util::DateTime subBegin(2023, 1, 1, 9, 0, 0);
    const util::Duration subLength("PT6H");
    eckit::LocalConfiguration configSubWindow;
    configSubWindow.set("begin", subBegin.toString());
    configSubWindow.set("length", subLength.toString());
    const util::TimeWindow subWindowBeginEnd(configSubWindow);
    // Expect this sub-window to have the same size as the first one.
    EXPECT(subWindowBeginEnd.start() == subWindow.start());
    EXPECT(subWindowBeginEnd.end() == subWindow.end());
    EXPECT(subWindowBeginEnd.length() == subWindow.length());
    EXPECT(subWindowBeginEnd.midpoint() == subWindow.midpoint());

    // Perform the same checks with values of `subBegin` and `subEnd` that exceed
    // the original window bounds.
    const util::DateTime subBeginLarge(2022, 12, 31, 0, 0, 0);
    const util::DateTime subEndLarge(2023, 1, 3, 0, 0, 0);
    const util::TimeWindow subWindowBeginEndLarge =
      timeWindow.createSubWindow(subBeginLarge, subEndLarge);
    EXPECT(subWindowBeginEndLarge.start() == subWindowLarge.start());
    EXPECT(subWindowBeginEndLarge.end() == subWindowLarge.end());
    EXPECT(subWindowBeginEndLarge.length() == subWindowLarge.length());
    EXPECT(subWindowBeginEndLarge.midpoint() == subWindowLarge.midpoint());
  }

  void testMask() {
    // Set up time window.
    const util::DateTime windowStart(2023, 1, 1, 0, 0, 0);
    const util::Duration windowLength("P1D");
    eckit::LocalConfiguration config;
    config.set("begin", windowStart.toString());
    config.set("length", windowLength.toString());
    const util::TimeWindow timeWindow(config);

    // Observed DateTimes.
    const std::vector<util::DateTime> obsTimes
      ({util::DateTime(2022, 12, 31, 12, 0, 0),
        util::DateTime(2023,  1,  1,  0, 0, 0),
        util::DateTime(2023,  1,  1, 12, 0, 0),
        util::DateTime(2023,  1,  2,  0, 0, 0),
        util::DateTime(2023,  1,  2, 12, 0, 0)});

    // Expected mask.
    // By default the TimeWindow object has an inclusive upper bound.
    const std::vector<bool> expectedMask({false, false, true, true, false});
    EXPECT(timeWindow.createTimeMask(obsTimes) == expectedMask);

    // Configure the TimeWindow object to use an inclusive lower bound and repeat the test.
    config.set("bound to include", "begin");
    const util::TimeWindow timeWindowInclusiveLower(config);
    const std::vector<bool> expectedMaskInclusiveLower
      ({false, true, true, false, false});
    EXPECT(timeWindowInclusiveLower.createTimeMask(obsTimes) ==
           expectedMaskInclusiveLower);

    // Perform the same tests after converting the DateTimes to a number of seconds
    // relative to an epoch.
    const util::DateTime epoch(2020, 1, 1, 0, 0, 0);
    std::vector<int64_t> obsEpochTimes;
    for (const auto & dt : obsTimes) {
      obsEpochTimes.push_back((dt - epoch).toSeconds());
    }
    // Check that an incorrect result is obtained if the epoch is not set.
    EXPECT_NOT(timeWindow.createTimeMask(obsEpochTimes) == expectedMask);
    EXPECT_NOT(timeWindowInclusiveLower.createTimeMask(obsEpochTimes) ==
               expectedMaskInclusiveLower);
    // Set the epoch correctly and repeat the test.
    timeWindow.setEpoch(epoch);
    EXPECT(timeWindow.createTimeMask(obsEpochTimes) == expectedMask);
    timeWindowInclusiveLower.setEpoch(epoch);
    EXPECT(timeWindowInclusiveLower.createTimeMask(obsEpochTimes) ==
           expectedMaskInclusiveLower);
  }

  CASE("util/TimeWindow/constructor") {
    testConstructor();
  }

  CASE("util/TimeWindow/accessors") {
    testAccessors();
  }

  CASE("util/TimeWindow/subwindow") {
    testSubWindow();
  }

  CASE("util/TimeWindow/mask") {
    testMask();
  }

class TimeWindow : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::TimeWindow";}
  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_TIMEWINDOW_H_
