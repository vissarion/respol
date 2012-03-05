# gnuplot script

set xlabel "Number of sequential convex hulls" font "Helvetica 30"
set ylabel "time"
#set yrange [0:2000]
set key right bottom
set logscale y

set terminal postscript clip lw 1 "Helvetica" 22
set output "hash-nohash-gfan.ps"

# TTR=traversing tropical resultant
# NFSI=normal fan from stable intersection

plot "results_gfan_d4/out_hash.txt" using 6:9 \
  smooth unique with linespoints title "hash", \
  "results_gfan_d4/out_nohash.txt" using 6:9 \
  smooth unique with linespoints title "no hash";
