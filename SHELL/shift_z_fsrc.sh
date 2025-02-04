#!/bin/sh -f

function zeropad {
expo=$((10 ** $2))
[ $1 -gt $expo ] &formatd=$(($1 + $expo))
echo ${formatd:1}
}

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
echo $nom
echo $num
NUMBER=`zeropad $num 3`
echo $NUMBER

cp $file fsrc.two
cp DATA/fshots fsrc.one
../BIN_LINUX/shift_fsrc
gnuplot <<EOF >$1/shift_z_fsrc.$NUMBER.gif
set term gif
set yrange [0:1.5]
plot 'fsrc.diz'
quit
EOF

done

cd $1

convert -delay 5 -loop 0 *.gif ../fsrc_z_anim.gif
#ffmpeg -i anim%05d.png anim.avi  
#animate  fsrc_anim.gif pour voir 
# ou bien firefox

cd ..




