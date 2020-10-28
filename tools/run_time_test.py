#!/usr/bin/env python3

"""
Comparison of two CTestCostData.txt based on
configurations in the given yaml file.

run_time_test.py ref_costdata.txt run_costdata.txt test_time.yaml

"""

import yaml
import sys

refFile = str(sys.argv[1])
runFile = str(sys.argv[2])
yamlFile = str(sys.argv[3])

with open(refFile) as f:
  ref_file = f.read()

with open(runFile) as f:
  run_file = f.read()

ref_file_dict = {}
for item in ref_file.split("\n"):
  if len(item.split(' ')) == 3:
    key = item.split(' ')[0]
    val = float(item.split(' ')[2])
    ref_file_dict[key] = val

run_file_dict = {}
for item in run_file.split("\n"):
  if len(item.split(' ')) == 3:
    key = item.split(' ')[0]
    val = float(item.split(' ')[2])
    run_file_dict[key] = val

with open(yamlFile) as f:
  #config_data = yaml.load(f, Loader=yaml.FullLoader)
  config_data = yaml.load(f)

status = "pass"

for tests in config_data:
  TestName = tests['TestName']
  Tol = tests['Tolerance']
  time_ref = ref_file_dict.get(TestName)
  time_run = run_file_dict.get(TestName)

  if all(t is not None for t in [time_ref, time_run]):
    rdiff = 100*(abs(float(time_ref) - float(time_run))/float(time_ref))

    if (time_ref < time_run  and rdiff > Tol):
      status = "fail"
      print("Test "+TestName+" failed")
      print("time_run = "+str(round(time_run,4)))
      print("time_ref = "+str(round(time_ref,4)))
      print("rdiff (%) = "+str(round(rdiff,2)))
    else:
      if status == "pass":
        status = "pass"
      print("Test "+TestName+" passed")
    print("***************************")
  else:
    print("time for "+TestName+" is missing")

if status == "fail":
  sys.exit(1)

