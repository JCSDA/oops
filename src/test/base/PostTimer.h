/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_BASE_POSTTIMER_H_
#define TEST_BASE_POSTTIMER_H_

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/base/PostTimer.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"

namespace test {

// -----------------------------------------------------------------------------
void testDefaultCtor() {
  const util::DateTime winbgn("2020-01-01T00:00:00Z");
  const util::DateTime winend("2020-01-02T00:00:00Z");
  const util::Duration step("PT1H");

  oops::PostTimer timer;
  timer.initialize(winbgn, winend);
  // any time between winbgn & winend should be fine
  EXPECT(timer.itIsTime(winbgn));
  EXPECT(timer.itIsTime(winbgn+step/2));
  EXPECT(timer.itIsTime(winbgn+step));
  EXPECT(timer.itIsTime(winend));
  EXPECT(!timer.itIsTime(winbgn-step));
  EXPECT(!timer.itIsTime(winend+step));
}

// -----------------------------------------------------------------------------
void testBgnEndCtor() {
  const util::DateTime winbgn("2020-01-01T00:00:00Z");
  const util::DateTime winend("2020-01-02T00:00:00Z");
  const util::Duration step("PT1H");

  // create a timer, make it run in winbgn->winend window
  oops::PostTimer timer1(winbgn, winend, step);
  // init with wider window (initial window should be used)
  timer1.initialize(winbgn-step, winend+step);
  // any time between winbgn & winend with step should be fine
  EXPECT(timer1.itIsTime(winbgn));
  EXPECT(!timer1.itIsTime(winbgn+step/2));
  EXPECT(timer1.itIsTime(winbgn+step));
  EXPECT(timer1.itIsTime(winend));
  EXPECT(!timer1.itIsTime(winbgn-step));
  EXPECT(!timer1.itIsTime(winend+step));

  // create a timer, make it run in winbgn-step->winend+step window
  oops::PostTimer timer2(winbgn-step, winend+step, step);
  // init with smaller window (initial window should be used)
  timer2.initialize(winbgn, winend);
  // any time between winbgn & winend with step should be fine
  EXPECT(timer2.itIsTime(winbgn));
  EXPECT(!timer2.itIsTime(winbgn+step/2));
  EXPECT(timer2.itIsTime(winbgn+step));
  EXPECT(timer2.itIsTime(winend));
  EXPECT(timer2.itIsTime(winbgn-step));
  EXPECT(timer2.itIsTime(winend+step));
}

// -----------------------------------------------------------------------------
void testConfCtor() {
  const util::DateTime winbgn("2020-01-01T00:00:00Z");
  const util::DateTime winend("2020-01-02T00:00:00Z");
  const util::Duration step("PT1H");

  eckit::LocalConfiguration test1;
  // create a timer with empty parameters (should be equiv to default ctor)
  oops::PostTimer timer1(oops::validateAndDeserialize<oops::PostTimerParameters>(test1));
  timer1.initialize(winbgn, winend);
  // any time between winbgn & winend should be fine
  EXPECT(timer1.itIsTime(winbgn));
  EXPECT(timer1.itIsTime(winbgn+step/2));
  EXPECT(timer1.itIsTime(winbgn+step));
  EXPECT(timer1.itIsTime(winend));
  EXPECT(!timer1.itIsTime(winbgn-step));
  EXPECT(!timer1.itIsTime(winend+step));

  const std::string freq = (step/2).toString();
  test1.set("frequency", freq);
  // create a timer, set frequency = step/2 through config
  oops::PostTimer timer2(oops::validateAndDeserialize<oops::PostTimerParameters>(test1));
  // init with winbgn->winend window
  timer2.initialize(winbgn, winend);
  // any time between winbgn & winend with step/2 should be fine
  EXPECT(timer2.itIsTime(winbgn));
  EXPECT(timer2.itIsTime(winbgn+step/2));
  EXPECT(!timer2.itIsTime(winbgn+step/4));
  EXPECT(timer2.itIsTime(winbgn+step));
  EXPECT(timer2.itIsTime(winend));
  EXPECT(!timer2.itIsTime(winbgn-step));
  EXPECT(!timer2.itIsTime(winend+step));

  eckit::LocalConfiguration test2;
  const std::string first = step.toString();
  test2.set("first", first);
  // create a timer, set first = step through config
  oops::PostTimer timer3(oops::validateAndDeserialize<oops::PostTimerParameters>(test2));
  // init with winbgn->winend window
  timer3.initialize(winbgn, winend);
  // any time between winbgn+step & winend should be fine
  EXPECT(!timer3.itIsTime(winbgn));
  EXPECT(!timer3.itIsTime(winbgn+step/2));
  EXPECT(timer3.itIsTime(winbgn+step));
  EXPECT(timer3.itIsTime(winbgn+step*2));
  EXPECT(timer3.itIsTime(winend));
  EXPECT(!timer3.itIsTime(winbgn-step));
  EXPECT(!timer3.itIsTime(winend+step));

  std::vector<std::string> vtimes{"2020-01-01T00:00:00Z", "2020-01-01T03:00:00Z",
                                  "2020-01-01T04:00:00Z", "2020-01-01T08:00:00Z",
                                  "2020-01-01T12:00:00Z"};
  eckit::LocalConfiguration test3;
  std::vector<std::string> times{"2020-01-01T00:00:00Z", "2020-01-01T04:00:00Z",
                                 "2020-01-01T08:00:00Z"};
  test3.set("times", times);
  std::vector<std::string> steps{"PT3H", "PT12H", "PT8H"};
  test3.set("steps", steps);
  // create a timer, set predefined steps
  oops::PostTimer timer4(oops::validateAndDeserialize<oops::PostTimerParameters>(test3));
  // init with winbgn->winend window
  timer4.initialize(winbgn, winend);
  // only times specified above should work
  EXPECT(timer4.itIsTime(winbgn));
  EXPECT(!timer4.itIsTime(winbgn+step));
  EXPECT(!timer4.itIsTime(winend));
  EXPECT(!timer4.itIsTime(winbgn-step));
  EXPECT(!timer4.itIsTime(winend+step));
  for (auto time : vtimes) {
    EXPECT(timer4.itIsTime(util::DateTime(time)));
  }
}

class PostTimer : public oops::Test {
 private:
  std::string testid() const override {return "test::PostTimer";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("util/PostTimer/defaultCtor") {
                      testDefaultCtor();
                    });
    ts.emplace_back(CASE("util/PostTimer/bgnEndCtor") {
                      testBgnEndCtor();
                    });
    ts.emplace_back(CASE("util/PostTimer/confCtor") {
                      testConfCtor();
                    });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_BASE_POSTTIMER_H_

