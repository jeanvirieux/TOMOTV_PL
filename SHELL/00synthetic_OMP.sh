# ############################### REFERENCE MODEL
#  NUMBER OF ITERATIONS MAXIMAL
MAX=1
I=1

# ###############################
#  COPY DATA
/bin/cp DATA/modelP.syn modelP
/bin/cp DATA/fsrc fsrc
/bin/cp DATA/fsta fsta
/bin/cp DATA/fobs fobs
/bin/cp DATA/fwei fwei	

/bin/cp DATA/model.head model.head
/bin/cp DATA/inversion.head inversion.head
/bin/cp DATA/cylinder.head .

echo " starting the tomography "



while [ $I -le $MAX ]
do
echo "---------------------------------------------------------"
echo "INTERPOLATED MODEL                ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/model.tomo


echo "---------------------------------------------------------"
echo "RAY TRACING                       ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/raytrace.tomo_OMP

echo "---------------------------------------------------------"
echo "DERIVATE                          ITERATION NUMBER"": "${I}
echo " "
../BIN_LINUX/derive.tomo_HFM

mv fcal DATA/fobs


echo $I
I=$((I+1))
done

# ######################## END LOOP




