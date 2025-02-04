#!/bin/sh -f

gnuplot <<EOF >rms_reduction.pdf
   set term pdf
   plot 'rms.file'
   quit
EOF


