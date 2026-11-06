/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_BASE_GETVALUESUTILS_H_
#define TEST_BASE_GETVALUESUTILS_H_

#include <string>
#include <vector>

// SonarQube wants this to be a constexpr, but that would be incompatible with eckit.
#define ECKIT_TESTING_SELF_REGISTER_CASES 0  // NOSONAR

#include "eckit/testing/Test.h"

#include "oops/base/ObsBufferDescriptor.h"
#include "oops/base/TimeInterpolationAccumulator.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/missingValues.h"

namespace test {

// -----------------------------------------------------------------------------
// ObsBufferDescriptor
// -----------------------------------------------------------------------------

void testLayoutSingleSamplingMethod() {
  // Two tasks, one sampling method, 4 levels per obs.
  const oops::ObsBufferDescriptor layout({{3}, {5}});
  const std::vector<size_t> varsizes{4};

  EXPECT(layout.ntasks() == 2);
  EXPECT(layout.nsms() == 1);
  EXPECT(layout.obsCount(0, 0) == 3);
  EXPECT(layout.obsCount(1, 0) == 5);

  EXPECT(layout.bufferSize(0, varsizes) == 12);
  EXPECT(layout.bufferSize(1, varsizes) == 20);
}

// -----------------------------------------------------------------------------

void testLayoutSumsSamplingMethods() {
  // One task, three sampling methods with different obs counts and variable sizes.
  const oops::ObsBufferDescriptor layout({{2, 3, 1}});
  const std::vector<size_t> varsizes{5, 2, 7};

  EXPECT(layout.bufferSize(0, varsizes) == 23);   // 2 * 5 + 3 * 2 + 1 * 7
}

// -----------------------------------------------------------------------------

void testLayoutEmptyBlocksAndTasks() {
  // A task may hold no observations at all for some, or all, sampling methods.
  const oops::ObsBufferDescriptor layout({{0, 0}, {0, 4}});
  const std::vector<size_t> varsizes{3, 2};

  EXPECT(layout.bufferSize(0, varsizes) == 0);
  EXPECT(layout.bufferSize(1, varsizes) == 8);
}

// -----------------------------------------------------------------------------

void testLayoutSizesDependOnVariableSet() {
  // The same observation grouping is used with two variable sets: the GeoVaLs variables
  // on the nonlinear path and the linear variables on the TL/AD path. The layout must
  // give the sizes for whichever set it is handed.
  const oops::ObsBufferDescriptor layout({{2, 3}});
  const std::vector<size_t> varsizes{4, 1};
  const std::vector<size_t> linsizes{1, 1};

  EXPECT(layout.bufferSize(0, varsizes) == 11);   // 2 * 4 + 3 * 1
  EXPECT(layout.bufferSize(0, linsizes) == 5);    // 2 * 1 + 3 * 1
}

// -----------------------------------------------------------------------------
// TimeInterpolationAccumulator: weighting rule
// -----------------------------------------------------------------------------

void testAccumulatorWeights() {
  const util::DateTime stateTime("2020-01-01T12:00:00Z");
  const util::Duration halfWidth("PT1H");
  const oops::TimeInterpolationAccumulator accumulator(stateTime, halfWidth);

  EXPECT(accumulator.isStateTime(stateTime));
  EXPECT(accumulator.weight(stateTime) == 0.0);

  // An observation after the state time is the state's first contribution, and is
  // weighted by how far it still is from the end of the interval.
  const util::DateTime after = stateTime + util::Duration("PT15M");
  EXPECT(accumulator.isFirstContribution(after));
  EXPECT(!accumulator.isStateTime(after));
  EXPECT(oops::is_close_absolute(accumulator.weight(after), 0.75, 1.0e-12));

  // An observation before the state time receives the state's second contribution.
  const util::DateTime before = stateTime - util::Duration("PT15M");
  EXPECT(!accumulator.isFirstContribution(before));
  EXPECT(oops::is_close_absolute(accumulator.weight(before), 0.75, 1.0e-12));

  // An observation at the edge of the interval receives no weight at all.
  EXPECT(oops::is_close_absolute(accumulator.weight(stateTime + halfWidth), 0.0, 1.0e-12));
  EXPECT(oops::is_close_absolute(accumulator.weight(stateTime - halfWidth), 0.0, 1.0e-12));
}

// -----------------------------------------------------------------------------

void testAccumulatorWeightsSumToOne() {
  // The two states bracketing an observation must contribute weights summing to one,
  // otherwise time interpolation would not reproduce a constant field.
  const util::Duration halfWidth("PT1H");
  const util::DateTime firstState("2020-01-01T12:00:00Z");
  const util::DateTime secondState = firstState + halfWidth;
  const util::DateTime obsTime = firstState + util::Duration("PT20M");

  const oops::TimeInterpolationAccumulator first(firstState, halfWidth);
  const oops::TimeInterpolationAccumulator second(secondState, halfWidth);

  EXPECT(first.isFirstContribution(obsTime));
  EXPECT(!second.isFirstContribution(obsTime));
  EXPECT(oops::is_close_absolute(first.weight(obsTime) + second.weight(obsTime),
                                 1.0, 1.0e-12));
}

// -----------------------------------------------------------------------------
// TimeInterpolationAccumulator: buffer walk
// -----------------------------------------------------------------------------

void testAccumulateAtStateTimeAssigns() {
  const util::DateTime stateTime("2020-01-01T12:00:00Z");
  const oops::TimeInterpolationAccumulator accumulator(stateTime, util::Duration("PT1H"));

  const std::vector obsTimes{stateTime};
  const std::vector mask{true};
  const std::vector<size_t> varLevels{2};
  const std::vector rhs{3.0, 4.0};
  std::vector lhs{99.0, 99.0};

  const size_t consumed = accumulator.accumulate(obsTimes, mask, varLevels, rhs, lhs.data());

  EXPECT(consumed == rhs.size());
  // An observation at the state time takes the state's value directly, unweighted, and
  // overwrites whatever was in the buffer.
  EXPECT(lhs[0] == 3.0);
  EXPECT(lhs[1] == 4.0);
}

// -----------------------------------------------------------------------------

void testAccumulateTwoContributionsBlend() {
  const util::Duration halfWidth("PT1H");
  const util::DateTime firstState("2020-01-01T12:00:00Z");
  const util::DateTime secondState = firstState + halfWidth;
  const util::DateTime obsTime = firstState + util::Duration("PT15M");

  const std::vector obsTimes{obsTime};
  const std::vector mask{true};
  const std::vector<size_t> varLevels{1};
  std::vector lhs{99.0};

  // The first state assigns its weighted value, discarding the buffer's previous content.
  const std::vector rhsFirst{10.0};
  oops::TimeInterpolationAccumulator(firstState, halfWidth)
    .accumulate(obsTimes, mask, varLevels, rhsFirst, lhs.data());
  EXPECT(oops::is_close_absolute(lhs[0], 7.5, 1.0e-12));

  // The second state adds to it.
  const std::vector rhsSecond{20.0};
  oops::TimeInterpolationAccumulator(secondState, halfWidth)
    .accumulate(obsTimes, mask, varLevels, rhsSecond, lhs.data());
  EXPECT(oops::is_close_absolute(lhs[0], 7.5 + 5.0, 1.0e-12));
}

// -----------------------------------------------------------------------------

void testAccumulateMaskedObsUntouched() {
  // Observations outside the state's time sub-window are masked out, and must keep the
  // values an earlier state gave them.
  const util::DateTime stateTime("2020-01-01T12:00:00Z");
  const oops::TimeInterpolationAccumulator accumulator(stateTime, util::Duration("PT1H"));

  const std::vector obsTimes{stateTime, stateTime};
  const std::vector mask{false, true};
  const std::vector<size_t> varLevels{2};
  const std::vector rhs{1.0, 2.0, 3.0, 4.0};
  const double sentinel = -777.0;
  std::vector lhs(4, sentinel);

  accumulator.accumulate(obsTimes, mask, varLevels, rhs, lhs.data());

  EXPECT(lhs[0] == sentinel);
  EXPECT(lhs[1] == sentinel);
  EXPECT(lhs[2] == 3.0);
  EXPECT(lhs[3] == 4.0);
}

// -----------------------------------------------------------------------------

void testAccumulateWalksMultiVariableBlock() {
  // Within a block the ordering is variable (slowest) -> obs (middle) -> level (fastest),
  // and variables may have different numbers of levels.
  const util::DateTime stateTime("2020-01-01T12:00:00Z");
  const oops::TimeInterpolationAccumulator accumulator(stateTime, util::Duration("PT1H"));

  const std::vector obsTimes{stateTime, stateTime};
  const std::vector mask{true, false};
  const std::vector<size_t> varLevels{3, 1};
  // Variable 0: obs 0 = {1,2,3}, obs 1 = {4,5,6}. Variable 1: obs 0 = {7}, obs 1 = {8}.
  const std::vector rhs{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  const double sentinel = -777.0;
  std::vector lhs(rhs.size(), sentinel);

  const size_t consumed = accumulator.accumulate(obsTimes, mask, varLevels, rhs, lhs.data());

  EXPECT(consumed == rhs.size());
  // Obs 0 is updated in both variables' sub-blocks...
  EXPECT(lhs[0] == 1.0);
  EXPECT(lhs[1] == 2.0);
  EXPECT(lhs[2] == 3.0);
  EXPECT(lhs[6] == 7.0);
  // ...and masked obs 1 is untouched in both.
  EXPECT(lhs[3] == sentinel);
  EXPECT(lhs[4] == sentinel);
  EXPECT(lhs[5] == sentinel);
  EXPECT(lhs[7] == sentinel);
}

// -----------------------------------------------------------------------------

void testAccumulateMissingValuesPropagate() {
  const double missing = util::missingValue<double>();
  const util::Duration halfWidth("PT1H");
  const util::DateTime firstState("2020-01-01T12:00:00Z");
  const util::DateTime secondState = firstState + halfWidth;
  const util::DateTime obsTime = firstState + util::Duration("PT15M");

  const std::vector obsTimes{obsTime};
  const std::vector mask{true};
  const std::vector<size_t> varLevels{2};

  // A missing interpolation result makes the observation missing...
  std::vector lhs{0.0, 0.0};
  const std::vector rhsFirst{missing, 10.0};
  oops::TimeInterpolationAccumulator(firstState, halfWidth)
    .accumulate(obsTimes, mask, varLevels, rhsFirst, lhs.data());
  EXPECT(lhs[0] == missing);
  EXPECT(oops::is_close_absolute(lhs[1], 7.5, 1.0e-12));

  // ...and it stays missing once the second state contributes, rather than being blended
  // with a partial value.
  const std::vector rhsSecond{20.0, 20.0};
  oops::TimeInterpolationAccumulator(secondState, halfWidth)
    .accumulate(obsTimes, mask, varLevels, rhsSecond, lhs.data());
  EXPECT(lhs[0] == missing);
  EXPECT(oops::is_close_absolute(lhs[1], 12.5, 1.0e-12));
}

// -----------------------------------------------------------------------------

class GetValuesUtils : public oops::Test {
 public:
  using oops::Test::Test;

 private:
  std::string testid() const override {return "test::GetValuesUtils";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("base/ObsBufferDescriptor/singleSamplingMethod") {
                      testLayoutSingleSamplingMethod();
                    });
    ts.emplace_back(CASE("base/ObsBufferDescriptor/sumsSamplingMethods") {
                      testLayoutSumsSamplingMethods();
                    });
    ts.emplace_back(CASE("base/ObsBufferDescriptor/emptyBlocksAndTasks") {
                      testLayoutEmptyBlocksAndTasks();
                    });
    ts.emplace_back(CASE("base/ObsBufferDescriptor/sizesDependOnVariableSet") {
                      testLayoutSizesDependOnVariableSet();
                    });
    ts.emplace_back(CASE("base/TimeInterpolationAccumulator/weights") {
                      testAccumulatorWeights();
                    });
    ts.emplace_back(CASE("base/TimeInterpolationAccumulator/weightsSumToOne") {
                      testAccumulatorWeightsSumToOne();
                    });
    ts.emplace_back(CASE("base/TimeInterpolationAccumulator/atStateTimeAssigns") {
                      testAccumulateAtStateTimeAssigns();
                    });
    ts.emplace_back(CASE("base/TimeInterpolationAccumulator/twoContributionsBlend") {
                      testAccumulateTwoContributionsBlend();
                    });
    ts.emplace_back(CASE("base/TimeInterpolationAccumulator/maskedObsUntouched") {
                      testAccumulateMaskedObsUntouched();
                    });
    ts.emplace_back(CASE("base/TimeInterpolationAccumulator/walksMultiVariableBlock") {
                      testAccumulateWalksMultiVariableBlock();
                    });
    ts.emplace_back(CASE("base/TimeInterpolationAccumulator/missingValuesPropagate") {
                      testAccumulateMissingValuesPropagate();
                    });
  }

  void clear() const override {
    // oops::Test requires this, but these cases hold no state between runs.
  }
};

}  // namespace test

#endif  // TEST_BASE_GETVALUESUTILS_H_
