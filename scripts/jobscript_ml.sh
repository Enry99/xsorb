#!/bin/bash
#Example of a jobscript for the machine learning optimization

# ...scheduler directives...

# Here you need to load/activate the python virtual env
# or conda environment where you have installed xsorb, so that
# the xsorb-ml-opt executable is available.
# The environment must also contain the libraries for the
# machine learning model you are using.


# Run the xsorb-ml-opt python executable
xsorb-ml-opt $1 $2 $3 $4