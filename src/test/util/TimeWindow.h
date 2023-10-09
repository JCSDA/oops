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

  void testAccessors() {
    // Set up time window.
    const util::DateTime windowStart(2023, 1, 1, 0, 0, 0);
    const util::DateTime windowEnd(2023, 1, 2, 0, 0, 0);
    const util::TimeWindow timeWindow(windowStart, windowEnd);

    // Test routines that access information about the time window.
    const util::DateTime expectedWindowStart(windowStart);
    const util::DateTime expectedWindowEnd(windowEnd);
    const util::Duration expectedWindowLength(expectedWindowEnd - expectedWindowStart);
    const util::DateTime expectedWindowMidpoint
      (expectedWindowStart + expectedWindowLength / 2);

    EXPECT(timeWindow.start() == expectedWindowStart);
    EXPECT(timeWindow.end() == expectedWindowEnd);
    EXPECT(timeWindow.length() == expectedWindowLength);
    EXPECT(timeWindow.midpoint() == expectedWindowMidpoint);

    // Create a TimeWindow object with an inclusive lower bound and repeat the tests.
    // Expect no difference in the results.
    const util::TimeWindow timeWindowInclusiveLower(windowStart, windowEnd,
                                                    util::InclusiveWindowBound::LOWER);
    EXPECT(timeWindowInclusiveLower.start() == expectedWindowStart);
    EXPECT(timeWindowInclusiveLower.end() == expectedWindowEnd);
    EXPECT(timeWindowInclusiveLower.length() == expectedWindowLength);
    EXPECT(timeWindowInclusiveLower.midpoint() == expectedWindowMidpoint);
  }

  void testTimeSlot() {
    // Set up time window.
    const util::DateTime windowStart(2023, 1, 1, 0, 0, 0);
    const util::DateTime windowEnd(2023, 1, 2, 0, 0, 0);
    const util::TimeWindow timeWindow(windowStart, windowEnd);

    // Valid time (in practice obtained from a model state).
    const util::DateTime validTime(2023, 1, 1, 12, 0, 0);
    // Window around valid time in which observations are accepted.
    const util::Duration halfWidth("PT3H");
    // Create time slot.
    const util::TimeWindow timeSlot = timeWindow.createSubWindow(validTime, halfWidth);

    // Check time slot accessors.
    const util::DateTime expectedTimeSlotStart(validTime - halfWidth);
    const util::DateTime expectedTimeSlotEnd(validTime + halfWidth);
    const util::Duration expectedTimeSlotLength(expectedTimeSlotEnd - expectedTimeSlotStart);
    const util::DateTime expectedTimeSlotMidpoint
      (expectedTimeSlotStart + expectedTimeSlotLength / 2);

    EXPECT(timeSlot.start() == expectedTimeSlotStart);
    EXPECT(timeSlot.end() == expectedTimeSlotEnd);
    EXPECT(timeSlot.length() == expectedTimeSlotLength);
    EXPECT(timeSlot.midpoint() == expectedTimeSlotMidpoint);

    // Perform the same checks with a large half-width. In this case expect the
    // time slot to have the same bounds as the original window.
    const util::Duration halfWidthLarge("P1D");
    const util::TimeWindow timeSlotLarge = timeWindow.createSubWindow(validTime, halfWidthLarge);

    EXPECT(timeSlotLarge.start() == timeWindow.start());
    EXPECT(timeSlotLarge.end() == timeWindow.end());
    EXPECT(timeSlotLarge.length() == timeWindow.length());
    EXPECT(timeSlotLarge.midpoint() == timeWindow.midpoint());
  }

  void testMask() {
    // Set up time window.
    const util::DateTime windowStart(2023, 1, 1, 0, 0, 0);
    const util::DateTime windowEnd(2023, 1, 2, 0, 0, 0);
    const util::TimeWindow timeWindow(windowStart, windowEnd);

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

    // Create a TimeWindow object with an inclusive lower bound and repeat the test.
    const util::TimeWindow timeWindowInclusiveLower(windowStart, windowEnd,
                                                    util::InclusiveWindowBound::LOWER);
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

  CASE("util/TimeWindow/accessors") {
    testAccessors();
  }

  CASE("util/TimeWindow/timeslot") {
    testTimeSlot();
  }

  CASE("util/TimeWindow/mask") {
    testMask();
  }

class TimeWindow : public oops::Test {
 private:
  std::string testid() const override {return "test::TimeWindow";}
  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_TIMEWINDOW_H_
