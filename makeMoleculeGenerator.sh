#/bin/sh
rm src/*.o
g++ -c -O3 -std=gnu++14 -o src/independentFunctions.o src/independentFunctions.cpp
g++ -c -O3 -std=gnu++14 -o src/point.o  src/point.cpp
g++ -c -O3 -std=gnu++14 -o src/polyHelix.o src/polyHelix.cpp 
g++ -c -O3 -std=gnu++14 -o src/randomMolGen.o src/randomMolGen.cpp
g++ -c -O3 -std=gnu++14 -o src/ktlMoleculeRandom.o src/ktlMoleculeRandom.cpp 
g++ -c -O3 -std=gnu++14 -o src/mainGlobalOptMultimerRanStart.o src/mainGlobalOptMultimerRanStart.cpp 
g++ -O3-std=gnu++14 -o makeMolecule src/independentFunctions.o src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/mainGlobalOptMultimerRanStart.o 
