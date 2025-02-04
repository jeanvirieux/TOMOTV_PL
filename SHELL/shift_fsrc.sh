#!/bin/sh -f

function zeropad {
expo=$((10 ** $2))
[ $1 -gt $expo ] &formatd=$(($1 + $expo))
echo ${formatd:1}
}

echo $1 

if [ $# -eq 0 ]
then
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  echo ' a parameter is needed: add h for help'
  echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  exit
fi

if [ $1 = h ]
then 

echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo directory where figures put for animation
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

else

  if [ $# -eq 1 ]
  then

rm fsrc.asc
rm fsrc.one
rm fsrc.two
rm fsrc.dis
rm fsrc.dix
rm fsrc.diy
rm fsrc.diz
rm fsrc.dat

mkdir $1

for file in fsrc.*

do

nom=`echo ${file%.*}`
num=`echo ${file#*.}`
#echo $nom
#echo $num
NUMBER=`zeropad $num 3`
echo $NUMBER

cp DATA/fshots fsrc.one
cp $file fsrc.two
../../BIN_LINUX/shift_fsrc

gnuplot <<EOF >$1/shift_d_fsrc.$NUMBER.gif
set term gif
set yrange [0:1.5]
plot 'fsrc.dis'
quit
EOF
gnuplot <<EOF >$1/shift_x_fsrc.$NUMBER.gif
set term gif
set yrange [0:1.5]
plot 'fsrc.dix'
quit
EOF
gnuplot <<EOF >$1/shift_y_fsrc.$NUMBER.gif
set term gif
set yrange [0:1.5]
plot 'fsrc.diy'
quit
EOF

gnuplot <<EOF >$1/shift_z_fsrc.$NUMBER.gif
set term gif
set yrange [0:1.5]
plot 'fsrc.diz'
quit
EOF

done

cd $1

convert -delay 5 -loop 0 shift_d_fsrc.*gif ../fsrc_d_anim.gif
#ffmpeg -i anim%05d.png anim.avi  
#animate  fsrc_anim.gif pour voir 
# ou bien firefox

convert -delay 5 -loop 0 shift_x_fsrc.*gif ../fsrc_x_anim.gif
convert -delay 5 -loop 0 shift_y_fsrc.*gif ../fsrc_y_anim.gif
convert -delay 5 -loop 0 shift_z_fsrc.*gif ../fsrc_z_anim.gif

cd ..
   else

echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo a parameter is needed: directory where to put figures for animation
echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   fi

fi


