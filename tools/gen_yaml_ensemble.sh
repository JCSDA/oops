#!/bin/bash

# Number of members required
nmem=10

# Copy first member here
cat<<EOFBASE >mem_0001
  - state:
    - date: 2010-01-01T00:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1D
    - date: 2010-01-01T01:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT1H
    - date: 2010-01-01T02:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT2H
    - date: 2010-01-01T03:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT3H
    - date: 2010-01-01T04:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT4H
    - date: 2010-01-01T05:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT5H
    - date: 2010-01-01T06:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT6H
    - date: 2010-01-01T07:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT7H
    - date: 2010-01-01T08:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT8H
    - date: 2010-01-01T09:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT9H
    - date: 2010-01-01T10:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT10H
    - date: 2010-01-01T11:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT11H
    - date: 2010-01-01T12:00:00Z
      filename: Data/forecast.ens.1.2009-12-31T00:00:00Z.P1DT12H
EOFBASE

# Prepare files
for imem in $(seq 2 ${nmem}) ; do
  imem4=`printf "%04d" ${imem}`
  sed -e s/"\.1\."/".${imem}."/g mem_0001 > mem_${imem4}
done

# Concatenate files
cat mem_* > mem_all

# Print all members
cat mem_all

# Delete all members
rm -f mem_*
