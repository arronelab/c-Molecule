#include "diffEvolRandom.h"

/***********************************************

  argv[1] should be the PDB name, e.g. 3V03
  argv[2] should be either "all", use all units or something like 12 (use units 1and  2 out of n_)
  argv[3] rand if a random start, nRand if not (currently st in the code)
  argv[4] should be the lower k bound 
  argv[5] should be the upper k bound
 **********************************************/

bool checkTransition(double &chiSqVal,double &chiSqCurr,double &uniformProb,int index,int &maxSteps){
  double tempFrac;
  if(index < maxSteps/3){
    tempFrac= 1.0-double(index)/double(maxSteps/3);
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
  }
}

bool checkTransitionSimple(double &chiSqVal,double &chiSqCurr,double &uniformProb,int index,int &maxSteps){
  if(chiSqVal<chiSqCurr){
    return true;
  }else{
    return false;
  }
}


int main( int argc, const char* argv[] )
{

  
 /*************************************

  set up model and scattering parameters 
 
  *************************************/
  double lmin=4.0; // closest distance two non adjactent local (same secondary unit) moelcules can get
  double rmin=3.7;double rmax=3.9; // max and min calpha-Dists
  double closestApproachDist=3.9; // closest distance two non adjactent non local moelcules (different secondary unit) can get
  ktlMolecule ktm;
  int noMols = std::atoi(argv[3]);
  ktm.setParams(rmin,rmax,lmin);
  ktm.readInMolGen(argv[1]);
  char fieldloc[1000]={}; 
  strcpy(fieldloc,argv[2]);
  int ind =1;
  const char* index  = intToChar(ind);
  strcat(fieldloc,index);
  strcat(fieldloc,".dat");
  ktm.getRandomMoleculeAllowOverlap();
  ktm.writeMoleculeToFile(fieldloc);
  for(int i=2;i<=noMols;i++){
 	ktm.getRandomMoleculeAllowOverlapReset();
	strcpy(fieldloc,argv[2]);
	index  = intToChar(i);
    	strcat(fieldloc,index);	
  	strcat(fieldloc,".dat");
	ktm.writeMoleculeToFile(fieldloc);
  }
  //ktm.removeOverlap();
  //ktm.writeMoleculeToFile(argv[2]);
  return 0;
}
      
