#! /bin/bash
#
# Runs my ADI PDE solver and creates a video from the
# output generated using gnuplot and ffmpeg.
#
# Author: Adam M. Holmes > /dev/null 2>&1


TIMESTEPS=1000
LAST=$(( $TIMESTEPS - 1 ))


### RUN PROGRAM
echo "[1/4] Running PDE solver..."
./heatpde $TIMESTEPS


### GENERATE PLOTS
echo "[2/4] Generating plots from data..."
for i in `seq 0 $LAST`;
do
TIME=$(echo "$i * 0.0002" | bc -l)
gnuplot << EOF
set term png
set output "data/$i.png"
set xrange [0:2]
set yrange [0:2]
set zrange [-1:1]
set title "exp(-t) * sin(πx) * sin(πy)\ntime: 0$TIME"
set pm3d at s
unset surface
unset colorbox
set palette rgbformulae 33,13,10
set nokey
splot 'data/$i.dat' title 'test' with lines
EOF
done


### CREATE VIDEO
echo "[3/4] Creating video from plots..."
rm heatpde.mp4 > /dev/null 2>&1
ffmpeg -start_number 0 -i data/%d.png heatpde.mp4 > /dev/null 2>&1


### CLEAN UP
echo "[4/4] Cleaning up..."
rm -rf data/*