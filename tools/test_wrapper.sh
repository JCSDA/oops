#!/bin/bash

# run the executable
echo ""
echo "==============================================================================="
echo "Running test executable"
echo "==============================================================================="

exename=$1
yamlname=$2
comparename=$3
runfile=$4
reffile=$5
ftol=$6
idif=$7
mpicmd=$8

# remove existing run file
rm testoutput/${runfile}

# Run Test
cmd="${mpicmd} ${exename} ${yamlname} testoutput/${runfile}"
echo ${cmd}
eval ${cmd}

exit_code=$?
if test "${exit_code}" == "0"; then
    echo -e "Test run succeed"
else
    echo -e "Test run failed with error code: ${exit_code} \n"
    exit ${exit_code}
fi

# Run compare
echo ""
echo "==============================================================================="
echo "Running compare script"
echo "==============================================================================="

cmd="${comparename}  testoutput/${runfile} testoutput/${reffile} ${ftol} ${idif}"
echo ${cmd}
eval ${cmd}

exit_code=$?
if test "${exit_code}" == "0"; then
   echo -e "Test compare succeed"
else
    echo -e "Test compare failed with error code: ${exit_code}"
    exit ${exit_code}
fi

# Test passed!
echo -e "PASSED"
