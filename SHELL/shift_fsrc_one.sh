#!/bin/sh -f

if [ $1 = h ]
then 

echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo premier parametre : iteration de la source
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

else

  if [ $# -eq 1 ]
  then

  cp fsrc.$1 fsrc.two
  cp DATA/fshots fsrc.one
  ../BIN_LINUX/shift_fsrc
gnuplot <<EOF >shift_fsrc.$1.pdf
   set term pdf
   plot 'fsrc.dis'
   quit
EOF

  else

  echo il faut un parametre faire h pour une aide

  fi

fi
