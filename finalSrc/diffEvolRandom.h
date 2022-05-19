#ifndef DIFF_EVOL
#define DIFF_EVOL

#include "scatteringCalculationRandom.h"
#include <random>
#include <ctime>
#include "independentFunctions.cpp"

class nelderMead{
 public:
  nelderMead(int maxStepsIn,double fitTolIn);
  nelderMead(int maxStepsIn,double fitTolIn,std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > &contactPairListIN);
  bool checkTransition(double &chiSqVal,double &chiSqCurr,double &uniformProb,int index);
  void minimizeLocalOneBestNLowMem(int &j,std::vector<scatteringCalculation> &sc);
  void minimizeLocalOneBestNLowMemWithContact(int &j,std::vector<scatteringCalculation> &sc,int index);
  void minimizeSmallVaritionNLowMemWithContact(int &j,std::vector<scatteringCalculation> &sc,int noIndiciesToChange,double &variationSize); 
  void minimizeLennyJ(int &j,std::vector<scatteringCalculation> &sc,int index);
  void minimizeLocalOneBestSingle(scatteringCalculation &sc);
  void minimizeLocalOneBestSingleWithContact(scatteringCalculation &sc);
  void minimizeLocalOneWriteToFile(int &j,std::vector<scatteringCalculation> &sc,std::ofstream &file);
  void minimizeGlobalWriteFile(std::vector<scatteringCalculation> &sc,const char* fieldname,double &kmin,double &kmax);
  void sortAndMixPop(std::vector<scatteringCalculation> &sc);
  void sortAndMixPopPair(std::vector<scatteringCalculation> &sc);
  void minimizeGlobalBestN(std::vector<scatteringCalculation> &sc,int &switchInt,int &writePeriod,const char* fieldname,double &kmin,double &kmax);
  void minimizeLocalBestN(scatteringCalculation &sc,int &writePeriod,const char* fieldname,double &kmin,double &kmax,double &prevBestIn);
  void writeBestFitToFile(const char* fieldname,scatteringCalculation &sc,int j);
  void writeBestFitToFileWPred(const char* fieldname,scatteringCalculation &sc);
  void writeBestFitToFileLocal(const char* fieldname,scatteringCalculation &sc);
  void writeBestFitToFileLocalWPred(const char* fieldname,scatteringCalculation &sc);
  void writeFinalFitsToFileWPred(const char* fieldname,std::vector<scatteringCalculation> &scx,int &index);
  void writeFinalFitsToFileWSbondPred(const char* fieldname,std::vector<scatteringCalculation> &sc,int &index);
  void writeLowLennyJToFile(const char* fieldname,std::vector<scatteringCalculation> &sc,int &index);
  void writeLowLennyJSbondsToFile(const char* fieldname,std::vector<scatteringCalculation> &sc,int &index);
  void writeConfigToFile(std::ofstream &file,std::pair<double,double> &chsqpr,std::vector<std::pair<double,double> > &ktpr);
  void minimizeLennyJPrePrep(std::vector<scatteringCalculation> &sc, int &maxLenJSteps);
  void minimizeGlobalBestNWContactPredictions(std::vector<scatteringCalculation> &sc,int &switchInt,int &writePeriod,const char* fieldname,double &kmin,double &kmax,int noTrialsBeforeSwitchToLocal,double localVariationSize);
  void minimizeLocalBestNWContactPredictions(scatteringCalculation &sc,int &writePeriod,const char* fieldname,double &kmin,double &kmax,std::pair<double,double> &initialBest);
 private:
  int maxSteps;
  double fitTol;
  std::vector<double> bestChiSq;
  std::vector<double> bestLennyJ;
  std::vector<std::pair<double,double> > bestChiSqPair;
  double bestChiSqSingle;
  std::pair<double,double> bestChiSqSinglePr;
  int changedIndex,molSize,end;
  std::vector<std::tuple<std::pair<int,int>,std::pair<int,int>,double> > contactPairList;
};


#endif
