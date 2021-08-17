#!/usr/bin/env python3

# (C) Copyright 2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

"""
Comparison of two text files containing Test log
from OOPS runs.

Comparison of floats uses a maximum relative difference:
(float1 - float2)/float1 <= tolerance.

Comparison of integers uses a maximum diffrence:
integer1 - integer2 <= difference

Comparison of strings is exact

Failure results in a return code of 1.

Call as:
compare.py run_file ref_file float_tolerance integer_difference
"""

import re
import sys

# Method that searches for int and floats in two lines
# and compares them one by one
def line_diff(line1,line2,lnum,ftol,idif):

  #Split line by whitespace or '='
  sline1 = re.split('\s+|=', line1)
  sline2 = re.split('\s+|=', line2)

  lineerror = 0

  for n in range(len(sline1)):
    if sline1[n] != '':

      nnm1 = rennm.findall(sline1[n])
      flt1 = reflt.findall(sline1[n])
      int1 = reint.findall(sline1[n])
      dat1 = redat.findall(sline1[n])

      nnm2 = rennm.findall(sline2[n])
      flt2 = reflt.findall(sline2[n])
      int2 = reint.findall(sline2[n])
      dat2 = redat.findall(sline2[n])

      found = 0

      #Compare if non-numeric string
      if nnm1 and nnm2:
        found=found+1
        if (nnm1 != nnm2):
          lineerror=lineerror+1
          print("Non numeric string mismatch at line "+str(lnum)+": "\
                +nnm1[0]+" not equal to "+nnm2[0])

      #Compare if string/float
      if flt1 and flt2:
        found=found+1
        flt1a = reflte.findall(flt1[0])
        flt2a = reflte.findall(flt2[0])
        rdiff = abs(float(flt1a[0])-float(flt2a[0]))/(abs(float(flt1a[0]))+1.0e-6)
        if (not rdiff <= ftol):
          lineerror=lineerror+1
          print("Float mismatch at line "+str(lnum)+": "+\
                flt1a[0]," not equal to ",flt2a[0]," with max relative difference ", ftol,\
                " Actual relative difference = ",rdiff)

      #Compare if integer
      if int1 and int2:
        found=found+1
        int1a = reinte.findall(int1[0])
        int2a = reinte.findall(int2[0])
        adiff = abs(int(int1a[0])-int(int2a[0]))
        if (not adiff <= idif):
          lineerror=lineerror+1
          print("Integer mismatch at line "+str(lnum)+": ",\
                int1a[0]+" not equal to "+int2a[0]," with max difference ", idif,\
                ". Actual difference = ",adiff)

      #Compare if date
      if dat1 and dat2:
        found=found+1
        if (dat1 != dat2):
          lineerror=lineerror+1
          print("Data string mismatch at line "+str(lnum)+": "+\
                dat1[0]+" not equal to "+dat2[0])


      #Compare if one is float and other is integer
      if flt1 and int2:
        found=found+1
        flt1a = reflte.findall(flt1[0])
        int2a = reinte.findall(int2[0])
        rdiff = abs(float(flt1a[0])-float(int2a[0]))/(abs(float(flt1a[0]))+1.0e-6)
        if (not rdiff <= ftol):
          lineerror=lineerror+1
          print("Float mismatch at line "+str(lnum)+": "+\
                flt1a[0]," not equal to ",int2a[0]," with max relative difference ", ftol,\
                " Actual relative difference = ",rdiff)

      if int1 and flt2:
        found=found+1
        int1a = reinte.findall(int1[0])
        flt2a = reflte.findall(flt2[0])
        rdiff = abs(float(int1a[0])-float(flt2a[0]))/(abs(float(int1a[0]))+1.0e-6)
        if (not rdiff <= ftol):
          lineerror=lineerror+1
          print("Float mismatch at line "+str(lnum)+": "+\
                int1a[0]," not equal to ",flt2a[0]," with max relative difference ", ftol,\
                " Actual relative difference = ",rdiff)

      #Exit with error if check has failed
      if found == 0:
        if sline1[n] != sline2[n]:
          print("In looping through line elements did not match either non numeric, float, date or integer. Error at line "\
                +str(lnum),". Trying to compare \'"+sline1[n]+"\' and \'"+sline2[n]+"\'")
          exit(1)

      if found > 1:
        if sline1[n] != sline2[n]:
          print("In looping through line elements matched multipe of non numeric, float, date and integer. Error at line "\
                +str(lnum)+". Trying to compare \'"+sline1[n]+"\' and \'"+sline2[n]+"\'")
          exit(1)

  return lineerror


# Get file names and tolerance from arguments
# -------------------------------------------
file_run = open(sys.argv[1], "r")
file_ref = open(sys.argv[2], "r")
ftol = float(sys.argv[3])
idif = int(sys.argv[4])

# Write grep results to a new file
file_runref = open(str(sys.argv[1])+".ref","w")

# Read reference file
lines_ref = file_ref.readlines()

# Potential regular expressions found in OOPS Test output
# -------------------------------------------------------
# Regex: ABC[12][ABC] (combination of string and int, e.g. AMSUA-NOAA19
rennm = re.compile('(^\D+[\d]*[\D]*$)')

# Regex: #[ABC][-]12.34[e[+-]12][,] (combination of string and flaot, e.g. MAX=123.123e-07,
reflt = re.compile('(^[\D]*?[-]?\d+\.\d+(?:[e][+-]?[\d]+)?[\,]?$)')   #[ABC][-]12.34[e[+-]12][,]

# Regex: #[-]12345[,] (combination of integer and comma, e.g. 123,
reint = re.compile('(^[-]?\d+[\,]?$)')

# Regex: #YYYY-MM-DDTHH:MN:SSZ[:] (date with potential semi-colon)
redat = re.compile('(^\d{4}[-]\d{2}[-]\d{2}[T]\d{2}[:]\d{2}[:]\d{2}[Z][\:]?$)')

# Sub extractions to convert combination of string and number to just number
reflte = re.compile('([+-]?\d+\.\d+(?:[e][+-]?[\d]+)?)')      #Float extraction (MAX=123.123e-07, -> 123.123e-07)
reinte = re.compile('([+-]?\d+)')                         #Integer extraction (123, -> 123)


# Loop through the run file and search on test string
refline = 0
error = 0
testfound = False
for line_run in file_run:
   if re.search('Test     : ', line_run):

       #Strip line
       line_strip = line_run[line_run.find('Test     : '):]

       #Check strings and integers
       lineerror = line_diff(line_strip,lines_ref[refline],refline+1,ftol,idif)
       error=error+lineerror

       #Write ref file in case update needed
       file_runref.write(line_strip)

       #Tick the reference line
       refline = refline + 1

       testfound = True

if (len(lines_ref) != refline):
  print("Test failed. "+str(refline-1)+" matches in run, "+str(len(lines_ref))+" ref")
  sys.exit(1) #Return failure

# Close the new reference file
file_runref.close()

# Return status
if error > 0:
  sys.exit(1) #Return failure
if not testfound:
  print("Did not find any instances of \'Test     : \' in run file")
  sys.exit(1) #Return failure

# Otherwise return success
sys.exit(0)
