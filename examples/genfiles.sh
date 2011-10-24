#!/bin/bash
#for i in `seq 10 29`;
#do
	#for j in `seq 0 1`;
	#	do
			for r in `seq 20 20 500`;
			do
				file="f"$r
			  echo  $file
			  cat experiments/restricted_hi/f
			  echo $r > experiments/restricted_hi/d8/$file
			  cat experiments/restricted_hi/f >> experiments/restricted_hi/d8/$file
			done
	#done
#done
		
