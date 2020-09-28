#!/usr/bin/env python3

# (C) Copyright 2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

"""
Codebuild Runs this script to check for branches with
same names accross downstream repos. Then creates a
webhook on PR page to notify users if such a branch exists.

Call as:
update_webhook_branchname.py branchName commitId
"""

import os
import boto3
from github import Github
import sys

branchName = sys.argv[1]
commitId = sys.argv[2]


# check other repos
repoList = ["ufo", "saber", "ioda", "soca", "fv3-jedi", "mpas-jedi", "shallow-water"]

token = os.getenv('GIT_PASS', '...')
g = Github(token)
owner = "JCSDA"

oopsRepo = g.get_repo("jcsda/oops")

for repoName in repoList:
  branchExists = True
  repo = g.get_repo(f"{owner}/{repoName}")
  branchList = list(repo.get_branches())
  try:
    repo.get_branch(branchName)
  except Exception:
    pass
    branchExists = False

  if (branchExists):
    print(branchName + ' exists in '+repoName)
    stageDescription= branchName+" exists in "+ repoName
    commitStatus = oopsRepo.get_commit(sha=commitId).create_status(
      state="pending",
      description=stageDescription,
      context="Branch Check-"+repoName)
  else:
    print(branchName + ' does not exist in '+repoName)
