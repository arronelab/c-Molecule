#/bin/sh
rm finalSrc/*.o
g++ -c -O3 -std=gnu++14 -o finalSrc/independentFunctions.o finalSrc/independentFunctions.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/point.o  finalSrc/point.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/polyHelix.o finalSrc/polyHelix.cpp 
g++ -c -O3 -std=gnu++14 -o finalSrc/randomMolGen.o finalSrc/randomMolGen.cpp
g++ -c -O3 -std=gnu++14 -o finalSrc/ktlMoleculeRandom.o finalSrc/ktlMoleculeRandom.cpp 
g++ -c -O3 -std=gnu++14 -o finalSrc/mainGlobalOptMultimerRanStart.o finalSrc/mainGlobalOptMultimerRanStart.cpp 
g++ -O3-std=gnu++14 -o makeMolecule finalSrc/independentFunctions.o finalSrc/point.o finalSrc/polyHelix.o finalSrc/randomMolGen.o finalSrc/ktlMoleculeRandom.o finalSrc/mainGlobalOptMultimerRanStart.o 
