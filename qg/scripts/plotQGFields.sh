#!/bin/bash

# Get back datelist
eval "datelist=(${datelist})"
ndate=${#datelist[@]}

if [ "${difference}" = true ] ; then
   eval "datelist_ref=(${datelist_ref})"
   ndate_ref=${#datelist_ref[@]}
   if [ ${ndate_ref} /= ${ndate} ] ; then
      echo "Inconsistent number of dates"
      exit 1
   fi
fi

if [ "${generate_gif}" = true ] ; then
   # Prepare frame directory
   mkdir -p tmp
fi

# Loop over dates
for ((idate=0;i<${ndate};idate++)); do
   # Export date
   export date=${datelist[$idate]}
   if [ "${difference}" = true ] ; then
      export date_ref=${datelist_ref[$idate]}
   fi
   echo "Running NCL script for "${date}

   # Run NCL
   ncl plotQGFields.ncl

   if [ "${generate_gif}" = true ] ; then
      # Link frames
      ln -sf `pwd`/fig/${filebase}.${date}.png `pwd`/tmp/`printf "%05d\n" ${idate}`.png
   fi

   # Update
   let "i++"
done

if [ "${generate_gif}" = true ] ; then
   # Generate gif
   ffmpeg -y -framerate 2 -i tmp/%05d.png fig/${filebase}.gif

   # Remove frames
   rm -fr tmp
fi
