# gnuplot script

#set title "Convex hull of points in a cube"
set xlabel "Number of points" font "Helvetica 30"
set ylabel "time"
set key right bottom
set logscale y 
#set xrange [0:200]

# 2:3 is offline
# 2:4 is online

plot "../../cube_2_hash.txt" using 2:4 \
    title "hash" smooth unique with lines lt 1, \
    "../../cube_2_nohash.txt" using 2:4 \
    title "no hash" smooth unique with lines lt 2, \
    "../../cube_3_hash.txt" using 2:4 \
    notitle smooth unique with lines lt 1, \
    "../../cube_3_nohash.txt" using 2:4 \
    notitle smooth unique with lines lt 2, \
    "../../cube_4_hash.txt" using 2:4 \
    notitle smooth unique with lines lt 1, \
    "../../cube_4_nohash.txt" using 2:4 \
    notitle smooth unique with lines lt 2, \
    "../../cube_5_hash.txt" using 2:4 \
    notitle smooth unique with lines lt 1, \
    "../../cube_5_nohash.txt" using 2:4 \
    notitle smooth unique with lines lt 2, \
    "../../cube_6_hash.txt" using 2:4 \
    notitle smooth unique with lines lt 1, \
    "../../cube_6_nohash.txt" using 2:4 \
    notitle smooth unique with lines lt 2, \
    "../../cube_7_hash.txt" using 2:4 \
    notitle smooth unique with lines lt 1, \
    "../../cube_7_nohash.txt" using 2:4 \
    notitle smooth unique with lines lt 2;

set terminal postscript eps color lw 1 "Helvetica" 22
#terminal postscript eps color lw 15 "Helvetica" 20
set output "ch_online_cube.ps"
replot
