#!/bin/bash

# parameters of the program are:
# $1: dimension
# $2: minimum number of points
# $3: maximum number of points
# $4: step in the number of points
if [ $# -lt 4 ]; then
        echo usage: $0 dimension numberofpoints
        exit 1
fi

d=$1
n=$2
n_max=$3
n_step=$4

echo \# D N T_h comp hashed T_nh
while [ $n -lt $n_max ]; do
        ./compare.sh $d $n
        let "n += $n_step"
done

exit 0

