#!/bin/bash
#for i in `seq 10 29`;
#do
	#for j in `seq 0 1`;
	#	do
			for r in `seq 20 20 500`;
			do
				file="cayley4_generic_d8"
				folder="gd8"
			  echo  $file$r
			  cat experiments/restricted_hi/$file
			  echo $r > experiments/restricted_hi/$folder/$file$r
			  cat experiments/restricted_hi/$file >> experiments/restricted_hi/$folder/$file$r
			done
	#done
#done
		
