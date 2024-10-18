#!/bin/bash
#Example of a quantum espresso job script
#For vasp it is similar, but without -input $1 >> $2

# ...scheduler directives...
# ...module load commands...
# ...environment variable settings...

# Run the pw.x executable
# the important part is $1 >> $2 for the input files
###############################################################
srun (...options...) pw.x (...options...) -input $1 >> $2
###############################################################