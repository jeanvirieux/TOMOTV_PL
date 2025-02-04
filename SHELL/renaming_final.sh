#!/bin/bash -f
echo renaming last models



for d in run*/
do
    [ -L "${d%/}" ] && continue
        echo "$d"
	cd $d

	for fichier in `ls modelP.??`
	do
          echo $fichier
	done
          echo $fichier 'selected'
        cp $fichier modelP.final

	for fichier in `ls modelS.??`
	do
          echo $fichier
	done
          echo $fichier 'selected'
        cp $fichier modelS.final

	for fichier in `ls fdift.sel.??`
        do
          echo $fichier
        done
        echo $fichier 'selected'
        cp $fichier fdift.sel.final

	for fichier in `ls fresi.??`
        do
          echo $fichier
        done
        echo $fichier 'selected'
        cp $fichier fresi.final

	for fichier in `ls fresw.??`
        do
          echo $fichier
        done
        echo $fichier 'selected'
        cp $fichier fresw.final

	for fichier in `ls fsrc.??`
        do
          echo $fichier
        done
        echo $fichier 'selected'
        cp $fichier fsrc.final



	cd ..
done
