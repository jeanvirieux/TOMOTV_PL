for fichier in `ls *STA.N.00038.txt`
do
    head -n 1 $fichier
    tail -1 $fichier
    echo '=========='
done
