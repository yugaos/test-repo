#!/bin/bash
## Software modules
module reset system
module load OpenFOAM/6-foss-2019b
#module load OpenFOAM/v2006-foss-2020a
#module load OpenFOAM/v2206-foss-2022a
source $FOAM_BASH
#
wmake

