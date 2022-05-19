#ifndef SCAT_CALC
#define SCAT_CALC

#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "point.h"
#include "hydrationShellRandom.h"
#include <random>
#include <chrono>  // for high_resolution_clock

class scatteringCalculation{
 public:
  bool is_digits(std::string str);
  std::string removeSpaces(std::string input);
  scatteringCalculation(const char* molName,double RinInp,double RoutInp,double Rshell,int ntrivsIn,double helixRatioIn,int solventsPerLinkIn,const char* filename,std::vector<int> &molIndicies,int isRand,int isContact,double &mutualDistCutOffIn,double &kminIN,double &kmaxIN,double &rmin,double &rmax,double &lmin,int &isSeqOrPDB);
  scatteringCalculation(const char* molName,double RinInp,double RoutInp,double RShell,int ntrivsIn,double helixRatioIn,int solventsPerLinkIn,const char* filename,std::vector<int> &molIndicies,int isRand,int isContact,double &mutualDistCutOffIn,double &kminIN,double &kmaxIN,double &rmin,double &rmax,double &lmin,const char *lenJFileName,int &WSeq);
  scatteringCalculation(){};
  double gaussian(double &s,std::vector<double> &gaussianVals);
  double gaussianFixed(double &s);
  double gaussianHydFixed(double &s);
  double hydrationScatter(double &s);
  int getNoScatPoints();
  int getMolSize();
  int getNoSections();
  void smoothDataLoggedRT(std::vector<std::pair<double,double> > molphase);
  void writeSmoothedDataToFile(const char* fieldname,const char* fileName);
  void writeModelScattering(const char* fieldname,const char* filename);
  void writeMoleculeToFile(const char* filename);
  void writeHelixPtsToFile(const char* filename);
  void writeHydrationLayerToFile(const char* filename);
  void writeHydrationParametersToFile(const char* filename);
  void writeRawScatteringToFile(const char* filename);
  std::vector<std::pair<double,double> > scatteringGaussianExcVolOrig();
  void printOneSelfScatter(int index);
  double chiSquaredExVolOrig();
  double chiSquaredExVolOrigTest();
  std::pair<double,double> chiSquaredExVolOrigWithContact(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  double lennardJonesOrig(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  double chiSquaredExVolReset();
  std::pair<double,double> chiSquaredExVolResetWithContact(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  double getPairDistance(std::pair<int,int> &index1,std::pair<int,int> &index2);
  double lennardJones(double &idealDist,double &currDist,int noConPred,double &weightCoeff);
  void getContactLists(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList);
  double chiSquaredAlter(int &index);
  std::pair<double,int> chiSquaredCollection(int &index);
  double chiSquaredAlterFast(int &index);
  void getScatteringBinnedRTUpdate();
  std::pair<double,double> chiSquaredAlterWithContact(int &index,std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  std::pair<double,double> chiSquaredAlterSmallVariationWithContact(std::vector<int> &changeIndicies,double &variationSize,std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  double alterLennyJ(int &index,std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  double getLennyJ(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  void setGaussianVals(std::vector<double> &gaussianVals);
  void setHydrationGaussianVals(std::vector<double> &gaussianVals);
  void setPhases(std::vector<std::pair<double,double> > solphase,std::vector<std::pair<double,double> > solmolphase,std::vector<std::pair<double,double> > molphase);
  void getBackboneStats();
  void changeMoleculeSingle(int changeIndex);
  void changeMoleculeSingle(int &index);
  void changeMoleculeSet(std::vector<int> &indicies);
  void changeMoleculeSingleMulti(int &index,int secIn);
  void changeMoleculeSetMulti(std::vector<int>  &indicies,int secIn);
  void changeMoleculeMultiRotate(double &angle,point &k,int secIn,point &transVec);
  void replicateMolecule(int &noReplications);
  void changeMoleculeLocal(int &index,double variationSize);
  void changeMoleculeLocalSet(std::vector<int> &indicies,double variationSize);
  void changeMoleculeSingle(int &index,std::vector<std::vector<point> > &cdsIn);
  int getSubsecSize(int sec);
  double maxNeighbourDistSec(int &sec);
  std::vector<double> checkOverlapWithRad(double &wRad);
  void getHydrophobicResidues();
  void getCoiledCoilResidues();
  double getGlobalRadiusOfCurvature();
  double getGlobalRadiusOfCurvatureBetweenSec();
  double coiledCoilPotential();
  double coiledCoilPotentialBetween();
  double getFitQualSbonds(std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairList,double weightCoeff);
  double compareDistances(std::vector<std::vector<point> > &coords2);
  std::string getType(int &chainIndex,int &index);
  std::vector<std::vector<point> > getCoordinates();
  void getScatteringBinnedRT(double &kmin,double &kmax);
  void printOnePod(int index);
  std::vector<std::pair<double,double> > getSolScatter();
  std::vector<std::pair<double,double> > getMolScatter();
  std::vector<std::pair<double,double> > getSolMolScatter();
  std::vector<std::vector<std::pair<double,double> > > getSelfSolScatter();
  std::vector<std::vector<std::pair<double,double> > > getSelfMolScatter();
  std::vector<std::vector<std::pair<double,double> > > getSelfSolMolScatter();
  std::vector<std::vector<std::vector<std::pair<double,double> > > > getMutualSolScatter();
  std::vector<std::vector<std::vector<std::pair<double,double> > > > getMutualMolScatter();
  std::vector<std::vector<std::vector<std::pair<double,double> > > > getMutualSolMolScatter();
  std::vector<std::pair<double,double> > getKapTauVals();
private:
  std::vector<point> rawData;
  std::vector<point> smoothedData;
  std::vector<int> emptyBins;
  std::vector<std::pair<double,double> > phasesMol;
  std::vector<std::pair<double,double> > phasesSol;
  std::vector<std::pair<double,double> > phasesSolMol;
  double a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,c;
  double a1H,a2H,a3H,a4H,a5H,b1H,b2H,b3H,b4H,b5H,cH;
  int nk;
  double kmin,kmax;
  hydrationShellMinimal hsmin;
};

#endif
