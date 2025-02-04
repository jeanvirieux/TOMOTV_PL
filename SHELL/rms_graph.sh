#!/bin/bash -f

awk '{print NR  "   " $1}' rms.pond > tt.dat
python ~/PYTHON/plot_xy_batch.py tt.dat
mv XY_2D.pdf rms_pond.pdf

awk '{print NR  "   " $1}' rms.file > tt.dat
python ~/PYTHON/plot_xy_batch.py tt.dat
mv XY_2D.pdf rms_file.pdf
