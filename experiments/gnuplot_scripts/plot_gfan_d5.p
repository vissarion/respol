# gnuplot script

set xlabel "Number of points" font "Helvetica 30"
set ylabel "time
set xrange [*:25]
#set yrange [0:2000]
set key bottom right
set logscale y

set terminal postscript clip lw 1 "Helvetica" 22
set output "hash-nohash-gfan-d5.ps"

# TTR=traversing tropical resultant
# NFSI=normal fan from stable intersection

plot "../results_gfan_d5/out_hash_d5.txt" using 4:9 \
  smooth unique with linespoints title "Respol-hash", \
  "../results_gfan_d5/out_nohash_d5.txt" using 4:9 \
  smooth unique with linespoints title "Respol-no hash", \
  "../results_gfan_d5/out_gfan_d5_2.txt" using 1:2 \
  smooth unique with linespoints pt 4 lt 4 title "Gfan-NFSI"

