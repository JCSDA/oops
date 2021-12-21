#!/bin/bash

# (C) Copyright 2020-2021 UCAR
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set -eux

echo "This is QG showing off"

cd $WORKDIR

$MPICMD $JEDIEXEC $1 2> stderr.$$.log 1> stdout.$$.log

cat stdout.$$.log

cat stderr.$$.log

echo "This is QG showing off"

