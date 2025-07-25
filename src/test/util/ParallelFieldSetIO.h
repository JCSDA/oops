/*
 * (C) Crown copyright 2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_PARALLELFIELDSETIO_H_
#define TEST_UTIL_PARALLELFIELDSETIO_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/runs/Test.h"
#include "oops/util/ParallelFieldSetIO.h"

namespace test {

// -------------------------------------------------------------------------------------------------

class TestParallelFieldSetIOFixture : private boost::noncopyable {
 public:
    static const std::vector<eckit::LocalConfiguration>& configs();
    static const std::vector<atlas::FunctionSpace>& nativeFunctionSpaces();
    static const std::vector<util::ParallelFieldSetIO>& parallelFieldSetIOs();

    static void reset();

 private:
    static TestParallelFieldSetIOFixture& getInstance();

    std::vector<eckit::LocalConfiguration> setupConfigs() const;
    std::vector<atlas::FunctionSpace> setupNativeFunctionSpaces() const;
    std::vector<util::ParallelFieldSetIO> setupParallelFieldSetIOs() const;

    TestParallelFieldSetIOFixture();

    ~TestParallelFieldSetIOFixture();

    std::unique_ptr<const std::vector<eckit::LocalConfiguration>> configs_;
    std::unique_ptr<const std::vector<atlas::FunctionSpace>> nativeFunctionSpaces_;
    std::unique_ptr<const std::vector<util::ParallelFieldSetIO>> parallelFieldSetIOs_;
};

// -------------------------------------------------------------------------------------------------

void setupFieldSets(atlas::FieldSet&, atlas::FieldSet&, const atlas::FunctionSpace&,
                    const eckit::LocalConfiguration&);

// -------------------------------------------------------------------------------------------------

void test(const size_t);

// -------------------------------------------------------------------------------------------------

class TestParallelFieldSetIO : public oops::Test {
 public:
    using oops::Test::Test;
 private:
    std::string testid() const override;
    void register_tests() const override;
    void clear() const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_UTIL_PARALLELFIELDSETIO_H_
