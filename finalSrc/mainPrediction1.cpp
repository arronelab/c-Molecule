#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h> 


/***********************************************
  argv[1] sequence file
  argv[2] initial prediction coords  file (can be empty)
  argv[3] prediction file 
 **********************************************/

bool checkTransition(double &chiSqVal,double &chiSqCurr,double &uniformProb,int index,int &maxSteps){
  /*double tempFrac;
  if(index < maxSteps/2){
    tempFrac= 1.0-double(index)/double(maxSteps/2);
  }else{
    tempFrac=0.0000000001;
  }
  double annealProb;
  double tempVal;
  if(chiSqVal<chiSqCurr){
    annealProb=1.0;
  }else{
    annealProb = std::exp(-(chiSqVal-chiSqCurr)/tempFrac);
  }
  //std::cout<<index<<" "<<maxSteps<<" "<<chiSqVal<<" "<<chiSqCurr<<" "<<annealProb<<" "<<uniformProb<<"\n";
  if(annealProb>uniformProb){
    return true;
  }else{
    return false;
  }*/
  if(chiSqVal<chiSqCurr){
    return true;
  }else{
    return false;
    }
}

int main( int argc, const char* argv[] )
{  
 /*************************************

  set up model parameters 
 
  *************************************/
  
  double lmin=4.0; // closest distance two non adjactent local (same secondary unit) moelcules can get
  double rmin=3.7;double rmax=3.9; // max and min calpha-Dists
  double closestApproachDist=3.9; // closest distance two non adjactent non local moelcules (different secondary unit) can get

  /*************************************
  
   determine initial model 

   *************************************/

  ktlMolecule mol;
  if(strcmp(argv[2],"no_initial_prediction")==0){
    // read in from sequernce and sc structure pred
    mol.readInSequence(argv[1],rmin,rmax,lmin);
    // generate random start (with no overlap)
    mol.getRandomMolecule();
    mol.writeMoleculeToFile("testMol.dat");
    // identify hydrophobic residues
    mol.getHydrophobicResidues();
  }else{
    // read in from sequernce and sc structure pred
    mol.readInSequence(argv[1],rmin,rmax,lmin);
    // read in coordinates
    mol.readInCoordinates(argv[2]);
    //here would be the filling in missing section routine
    // identify hydrophobic residues
    mol.writeMoleculeToFile("testMol.dat");
    mol.getHydrophobicResidues();
  }
  bool doAll = false;
  
  mol.writeMoleculeToFile(argv[3]);