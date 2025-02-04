# ############################### REFERENCE MODEL
#  NUMBER OF ITERATIONS MAXIMAL
MAX=6
I=1

# ###############################
echo " continuing the tomography "



while [ $I -le $MAX ]
do
echo "---------------------------------------------------------"
echo "INTERPOLATED MODEL                ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/model.tomo


echo "---------------------------------------------------------"
echo "RAY TRACING                       ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/raytrace.tomo

echo "---------------------------------------------------------"
echo "DERIVATE                          ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/derive.tomo

cp fmat.x fmat.x_bak
cp fmat.ic fmat.ic_bak
cp fmat.id fmat.id_bak

echo "---------------------------------------------------------"
echo "COND.                             ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/precond.tomo

echo "------------------------- -------------------------------"
echo "INVERSION                         ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/inverse.tomo

echo "------------------------- -------------------------------"
echo "INVARIANCE                        ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/modelP2_5D

if [ $I -ge 1 ]
then
  /bin/cp modelP modelP.$I
  /bin/cp modelP.2d modelP_2d.$I
fi


echo $I
I=$((I+1))
done

# ######################## END LOOP

