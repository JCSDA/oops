/*
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_TEMPLATEDENSEMBLEAPPLICATION_H_
#define OOPS_RUNS_TEMPLATEDENSEMBLEAPPLICATION_H_

#include <string>
#include <vector>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"


namespace oops {


template <typename APP>
class TemplatedEnsembleApplication : public Application {
 public:
// -----------------------------------------------------------------------------

  explicit TemplatedEnsembleApplication(const eckit::mpi::Comm & comm = oops::mpi::world()) :
      Application(comm) {}

// -----------------------------------------------------------------------------

  virtual ~TemplatedEnsembleApplication() {}

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    // Define member applications from a template yaml and templating information

    const auto templatingInformation = fullConfig.getSubConfiguration("templating information");

    // Number of independent runs of the application
    const size_t nmembers = templatingInformation.getUnsigned("number of members");

    // Member index for the first member that is run in this executable; for example,
    // if `number of members` = 4 and `first member index` = 17, then the
    // ensemble members for this run have indices 17, 18, 19 and 20.
    const size_t firstMember = templatingInformation.getUnsigned("first member index", 1);

    // Maximum number of digits in the ensemble members' indices (the default choice of 3 digits
    // means that ensemble member indices cannot exceed 999).
    const size_t patternLength = templatingInformation.getUnsigned(
                                    "max number of digits in member index", 3);

    const size_t lastMember = firstMember - 1 + nmembers;
    const size_t ntasks = this->getComm().size();
    const size_t mytask = this->getComm().rank();
    const size_t tasks_per_member = ntasks / nmembers;
    const size_t mymember = firstMember + mytask / tasks_per_member;

    Log::info() << "Running members "
                << firstMember << " to " << lastMember << ", handled by " << ntasks
                << " total MPI tasks and " << tasks_per_member
                << " MPI tasks per member." << std::endl;

    ASSERT(ntasks%nmembers == 0);
    ASSERT(lastMember < std::pow(10.0, patternLength));

    //  Create the communicator for each member, named comm_member_{i}:
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = this->getComm().split(mymember, commName);

    //  Each member uses a different configuration:
    const auto memberTemplate = fullConfig.getSubConfiguration("template yaml");
    const auto memberConf = fillInTemplate(memberTemplate, templatingInformation, mymember);
    APP ensapp(commMember);
    return ensapp.execute(memberConf, validate);
    return 1;
  }

// -----------------------------------------------------------------------------

  void validateConfig(const eckit::Configuration & fullConfig) const override {
    if (!fullConfig.has("templating information")) {
        throw eckit::UserError("Missing required yaml entry `templating information`",
                               Here());
    }

    if (!fullConfig.has("template yaml")) {
        throw eckit::UserError("Missing required yaml entry `template yaml`",
                               Here());
    }

    const auto templatingInformation = fullConfig.getSubConfiguration("templating information");

    if (!templatingInformation.has("number of members")) {
        throw eckit::UserError("Missing required yaml parameter "
                               "`templating information: number of members",
                               Here());
    }

    if (!templatingInformation.has("pattern with zero padding")) {
        throw eckit::UserError("Missing required yaml parameter "
                               "`templating information: pattern with zero padding",
                               Here());
    }

    if (!templatingInformation.has("pattern without zero padding")) {
        throw eckit::UserError("Missing required yaml parameter "
                               "`templating information: pattern without zero padding",
                               Here());
    }

    // For ensemble applications we also need to validate individual yamls
    APP ensapp(oops::mpi::world());
    const auto memberTemplate = fullConfig.getSubConfiguration("template yaml");
    const size_t nMembers = templatingInformation.getUnsigned("number of members");

    for (size_t jj = 0; jj < nMembers; ++jj) {
      const auto memberConf = fillInTemplate(memberTemplate, templatingInformation, jj);
      ensapp.validateConfig(memberConf);
    }
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const override {
    return "oops::TemplatedEnsembleApplication<>";
  }

// -----------------------------------------------------------------------------

  eckit::LocalConfiguration fillInTemplate(
          const eckit::LocalConfiguration & memberTemplate,
          const eckit::LocalConfiguration & templatingInformation,
          const size_t member) const {
    // Template pattern in the member template yaml, to be substituted by the individual
    // ensemble members' indices left-padded with zeros (subject to max number of digits
    // in member index).
    const auto & patternWithPad = templatingInformation.getString(
                                    "pattern with zero padding");


    // Template pattern in the member template yaml, to be substituted by the individual
    // ensemble members' indices without zero-padding.
    const auto & patternNoPad = templatingInformation.getString(
                                    "pattern with zero padding");

    const auto patternLength = templatingInformation.getUnsigned(
                                    "max number of digits in member index", 3);

    eckit::LocalConfiguration memberConfig(memberTemplate);
    util::seekAndReplace(memberConfig, patternWithPad, member, patternLength);
    util::seekAndReplace(memberConfig, patternNoPad, std::to_string(member));
    return memberConfig;
  }

// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_TEMPLATEDENSEMBLEAPPLICATION_H_
