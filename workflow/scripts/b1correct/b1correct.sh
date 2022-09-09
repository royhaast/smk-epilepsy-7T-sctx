#!usr/bin/bash

: '
This script corrects the MP2RAGE data for B1+ inhomogeneities and creates 
a synthetic FLAWS image as in https://doi.org/10.1097/RLI.0000000000000718

Requirements:
- MATLAB
- CRMBM MP2RAGE processing script
'

script=$1
script_dir=`readlink -f $script`
script_dir=$(dirname "${script_dir}")
script=`basename $script .m`

out=$2

pushd $script_dir
    matlab -nosplash -nodisplay -nodesktop -r "${script}('${out}'); quit"
popd
