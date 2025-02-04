# ############################### REFERENCE MODEL
#  NUMBER OF ITERATIONS MAXIMAL
MAX=125
I=1

# ###############################
#  COPY DATA
/bin/cp DATA/modelP.ini modelP
/bin/cp DATA/modelS.ini modelS
/bin/cp DATA/fsrc fsrc
/bin/cp DATA/fsta fsta
/bin/cp DATA/fobs fobs
/bin/cp DATA/fwei fwei	

/bin/cp DATA/model.head model.head
/bin/cp DATA/inversion.head inversion.head

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

if [ $I -ge 1 ]
then
  /bin/cp modelP modelP.$I
  /bin/cp fsrc fsrc.$I
fi


echo $I
I=$((I+1))
done

# ######################## END LOOP




