#!/bin/sh

# the program receives two parameters, first is the dimension and second is
# the number of points
if [ $# -lt 2 ]; then
        echo usage: $0 dimension numberofpoints
        exit 1
fi

D=$1
N=$2

echo -n "$D $N "

# run hashed triangulation, output the time and store the stderr output in
# a file
rm -f out2
./triangulation_hash $D $N 2> out2
echo -n " "
cat out2 | ./filter_hash_output.pl
echo -n " "

# run non-hashed triangulation and output the time
./triangulation_no_hash $D $N
echo

rm -f out2

exit 0
