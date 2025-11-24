/*
 * (C) Crown copyright 2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "test/util/ParallelFieldSetIO.h"

#if !defined(AOCC) && !defined(NVHPC)
#include <filesystem>
#endif
#include <algorithm>

#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/util/function/VortexRollup.h"
#include "oops/runs/Run.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"

namespace test {

using Test_ = TestParallelFieldSetIOFixture;

// -------------------------------------------------------------------------------------------------

const std::vector<eckit::LocalConfiguration>& Test_::configs() {
    return *getInstance().configs_;
}
const std::vector<atlas::FunctionSpace>& Test_::nativeFunctionSpaces() {
    return *getInstance().nativeFunctionSpaces_;
}
const std::vector<util::ParallelFieldSetIO>& Test_::parallelFieldSetIOs() {
    return *getInstance().parallelFieldSetIOs_;
}

void Test_::reset() {
    getInstance().configs_.reset();
    getInstance().nativeFunctionSpaces_.reset();
    getInstance().parallelFieldSetIOs_.reset();
}

Test_& Test_::getInstance() {
    static Test_ theTest_;
    return theTest_;
}

std::vector<eckit::LocalConfiguration> Test_::setupConfigs() const {
    return TestEnvironment::config().getSubConfigurations("test configurations");
}

std::vector<atlas::FunctionSpace> Test_::setupNativeFunctionSpaces() const {
    std::vector<atlas::FunctionSpace> nativeFunctionSpaces;
    nativeFunctionSpaces.reserve(configs_->size());
    for (const auto& config : *configs_) {
        atlas::Grid grid{};
        atlas::grid::Partitioner partitioner{};
        atlas::Mesh mesh{};
        nativeFunctionSpaces.emplace_back(atlas::FunctionSpace{});
        atlas::FieldSet fieldSet{};
        util::setupFunctionSpace(atlas::mpi::comm(), config, grid, partitioner, mesh,
                                 nativeFunctionSpaces.back(), fieldSet);
    }
    return nativeFunctionSpaces;
}

std::vector<util::ParallelFieldSetIO> Test_::setupParallelFieldSetIOs() const {
    std::vector<util::ParallelFieldSetIO> parallelFieldSetIOs;
    parallelFieldSetIOs.reserve(configs_->size());
    for (size_t c = 0; c < configs_->size(); c++) {
        const std::string gridName = (*configs_)[c].getSubConfiguration("grid").getString("name");
        parallelFieldSetIOs.emplace_back((*nativeFunctionSpaces_)[c], gridName);
    }
    return parallelFieldSetIOs;
}

Test_::TestParallelFieldSetIOFixture()
    : configs_(std::make_unique<std::vector<eckit::LocalConfiguration>>(setupConfigs())),
    nativeFunctionSpaces_(
        std::make_unique<std::vector<atlas::FunctionSpace>>(setupNativeFunctionSpaces())),
    parallelFieldSetIOs_(
        std::make_unique<std::vector<util::ParallelFieldSetIO>>(setupParallelFieldSetIOs())) {}

Test_::~TestParallelFieldSetIOFixture() {}

// -------------------------------------------------------------------------------------------------

void setupFieldSets(atlas::FieldSet& source, atlas::FieldSet& target,
                    const atlas::FunctionSpace& fSpace, const eckit::LocalConfiguration& config) {
    const auto nFields = config.getUnsigned("nfields");
    std::vector<std::string> fNames(nFields, "field_");
    auto f = 0;
    std::transform(fNames.begin(), fNames.end(), fNames.begin(),
                   [&f](auto s) { return s + std::to_string(f++); });

    const std::vector<size_t> r2Sizes = config.getUnsignedVector("rank 2 sizes");
    const std::vector<size_t> r3Sizes = config.getUnsignedVector("rank 3 sizes");
    ASSERT_MSG(r2Sizes.size() == nFields,
               "test::setupFieldSets: no. rank 2 sizes does not match nfields");
    ASSERT_MSG(r3Sizes.size() == nFields,
               "test::setupFieldSets: no. rank 3 sizes does not match nfields");

    source = util::createSmoothFieldSet(atlas::mpi::comm(), fSpace, r2Sizes, r3Sizes, fNames);
    source.haloExchange();
    target = util::createFieldSet(fSpace, r2Sizes, r3Sizes, fNames);
}

// -------------------------------------------------------------------------------------------------

void test(const size_t configNumber) {
    const eckit::LocalConfiguration& config = Test_::configs()[configNumber];
    const atlas::FunctionSpace& nativeFunctionSpace = Test_::nativeFunctionSpaces()[configNumber];
    const util::ParallelFieldSetIO& parallelFieldSetIO = Test_::parallelFieldSetIOs()[configNumber];

    atlas::FieldSet source, target;
    setupFieldSets(source, target, nativeFunctionSpace, config);

    // Test symmetry of write and read methods on constant number of PEs
    const auto ncfilepath = TestEnvironment::config().getString("datadir") + "/" +
                            "test_util_parallelfieldsetio_" +
                            std::to_string(atlas::mpi::comm().size()) + "PE-" +
                            "config" + std::to_string(configNumber) + ".nc";

    if (source.size() == 1) {
        // If nfields == 1, check Field API
        parallelFieldSetIO.write(source[0], ncfilepath);
        parallelFieldSetIO.read(target[0], ncfilepath);
    } else {
        // Otherwise, check FieldSet API
        parallelFieldSetIO.write(source, ncfilepath);
        parallelFieldSetIO.read(target, ncfilepath);
    }
    EXPECT(util::compareFieldSets(source, target, 1.0e-16));

    if (atlas::mpi::comm().size() != 1) {
        // Test invariance of file to domain decomposition
        const auto ncfilepathSerial = TestEnvironment::config().getString("datadir") + "/" +
                                      "test_util_parallelfieldsetio_1PE-" +
                                      "config" + std::to_string(configNumber) + ".nc";
#if !defined(AOCC) && !defined(NVHPC)
        if (std::filesystem::exists(ncfilepathSerial)) {
#else
        if (access(ncfilepathSerial.c_str(), 0) == 0) {
#endif
            parallelFieldSetIO.read(target, ncfilepath);
            EXPECT(util::compareFieldSets(source, target, 1.0e-16));
        } else {
            oops::Log::warning() <<
                "util::ParallelFieldSetIO domain decomposition invariance test skipped as files "
                "produced with 1PE do not exist. Are test dependencies set up correctly?"
                                 << std::endl;
        }
    }
}

// -------------------------------------------------------------------------------------------------

CASE("test") {
    for (size_t c = 0; c < Test_::configs().size(); c++) {
        SECTION("config" + std::to_string(c)) {
            test(c);
        }
    }
}

// -------------------------------------------------------------------------------------------------

std::string TestParallelFieldSetIO::testid() const {
    return "oops::test::TestParallelFieldSetIO";
}
void TestParallelFieldSetIO::register_tests() const {}
void TestParallelFieldSetIO::clear() const {}

// -------------------------------------------------------------------------------------------------

}  // namespace test

int main(int argc,  char ** argv) {
    oops::Run run(argc, argv);
    test::TestParallelFieldSetIO tests;
    return run.execute(tests);
}
