#!/bin/bash
# ----------------------------------------------------------------------------------
# User-defined parameters
# ----------------------------------------------------------------------------------
# Data directory
datadir="${HOME}/build/ufo-bundle/oops/qg/test/testoutput"

# File bases
declare -a filebases=( "3densvar"
                       "3dvar"
                       "3dvar_change_var"
                       "3dvar_hybrid"
                       "4densvar"
                       "4dvar_dripcg"
                       "4dvar_drpcg_lmp"
                       "4dvar_drpfom"
                       "4dvar_drplanczos"
                       "4dvar_drplanczos_hybrid"
                       "4dvar_forcing"
                       "4dvar_ipcg"
                       "4dvar_obs_biased"
                       "4dvar_rpcg"
                       "4dvar_saddlepoint"
                     )
# ----------------------------------------------------------------------------------
# End of user-defined parameters
# ----------------------------------------------------------------------------------
nfilebases=${#filebases[@]}

for ((ifilebases=0;ifilebases<${nfilebases};ifilebases++)); do
   # Export parameters
   export datadir=${datadir}
   export filebase=${filebases[$ifilebases]}
   echo "Running NCL script for "${filebase}

   # Get cost function values
   mkdir -p tmp
   grep " CostJb   : Nonlinear Jb" ${datadir}/${filebase}.test.log.out | awk '{print $8}' > tmp/Jb
   grep " CostJo   : Nonlinear Jo" ${datadir}/${filebase}.test.log.out | awk '{print $8}' | sed 's/,//' > tmp/Jo
   grep " CostJo   : Nonlinear Jo" ${datadir}/${filebase}.test.log.out | awk '{print $14}' | sed 's/,//' > tmp/Jo_n
   grep " CostJcDFI: Nonlinear Jc" ${datadir}/${filebase}.test.log.out | awk '{print $7}' > tmp/Jc
   sed -n '/^Jo Departures/,/^End Jo Departures/p;/^End Jo Departures/q' ${datadir}/${filebase}.test.log.out | sed -e '1d' -e '$ d' | awk '{print $1}' > tmp/obstypes
   grep "max iter = " ${datadir}/${filebase}.test.log.out | awk '{print $5}' | sed 's/,//' > tmp/ninner
   grep "Quadratic cost function: Jb" ${datadir}/${filebase}.test.log.out | sed 's/^.*\( = .*\).*$/\1/' | awk '{print $2}' > tmp/Jb_quad
   grep "Quadratic cost function: JoJc" ${datadir}/${filebase}.test.log.out | sed 's/^.*\( = .*\).*$/\1/' | awk '{print $2}' > tmp/JoJc_quad

   # Run NCL
   ncl plotQGCost.ncl

   # Clean
   rm -fr tmp
done
