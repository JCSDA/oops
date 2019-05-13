#!/bin/bash
# ----------------------------------------------------------------------------------
# User-defined parameters
# ----------------------------------------------------------------------------------
# Data directory
datadir="/Users/cpaessvisitor/jedivm/vagrant_data/build/ufo-bundle/bin"

# File base
filebase="example.an"

# Date list
declare -a datelist=("2010-01-01T12:00:00Z"
                    )

# Reference file base
filebase_ref="example.fc"

# Reference date list
declare -a datelist_ref=("2009-12-31T00:00:00Z.P1DT12H"
                    )

# Gif option
generate_gif=false
# ----------------------------------------------------------------------------------
# End of user-defined parameters
# ----------------------------------------------------------------------------------
# Export parameters
export datadir=${datadir}
export filebase=${filebase}
export difference=true
export filebase_ref=${filebase_ref}
export generate_gif=${generate_gif}

# Run plotQGFields.sh
datelist="${datelist[@]@Q}" datelist_ref="${datelist_ref[@]@Q}" ./plotQGFields.sh
