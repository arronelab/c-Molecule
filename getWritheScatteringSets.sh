#!/bin/bash
ScatterFile=newFitData/$1/$1Saxs.dat
fileLocs=newFitData/$1/
initialCoordsFile=newFitData/$1/coordinates.dat
#initialCoordsFile=calphaData/1LYZ/1LYZch_1nameList.dat
noStructures=1
pairedPredictions=none
fixedsections=newFitData/$1/varyingSectionSecondary.dat 
crystalSymmetry=none
withinMonomerHydroCover=none
betweenMonomerHydroCover=none
kmin=0.01;
kmax=0.5;
maxNoFitSteps=2

mkdir newFitData/$1/$2


for i in {1..1}
do
    writheScatterTest $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps newFitData/$1/$2/mol$i newFitData/$1/$2 newFitData/$1/mixtureFile.dat
done
 
