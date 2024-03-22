#ifndef trekG4Mass2_h
#define trekG4Mass2_h 1
// Template class to calculate beta2 and mass2

#include <TMath.h>
#include <TRandom3.h>
#include <random>
#include "Pythia8/Pythia.h"
using namespace Pythia8;

template<typename T> 
class disvars{
public:
  disvars();
  disvars(Vec4 pProton,Vec4 peIn,Vec4 peOut){
  pProtonV4=pProton;
  peInV4=peIn;
  peOutV4=peOut;
  pPhoton=peIn-peOut;
  };
  inline void setJet4Vec(Vec4 jet4V){jet4vec=jet4V;};
  bool isBHadron(int id);
  double Q2();
  double W2();
  double xBjorken();
  double yLine();
  double eprimeTheta();
  double eprimeEnergy();
  double eprimepT();
  double myRapidity(Vec4 p);
  double jetTheta();
  double jetPhi();
  double jetEnergy();
  double jetpT();
  double theta_jet, E_jet, phi_jet;
protected:
  const double c=2.99792458e8;
private:
  // four-vector of proton, incoming and scattered e-
  Vec4 pPhoton,pProtonV4, peInV4, peOutV4, jet4vec;
  double tof1, tof2, plen, P, dt;
};
template<typename T> 
double disvars<T>::Q2(){
  return -1*pPhoton.m2Calc();
}
template<typename T> 
double disvars<T>::W2(){
  return (pProtonV4+pPhoton).m2Calc();
}
template<typename T> 
double disvars<T>::eprimeTheta(){
  return 180/TMath::Pi()*(peOutV4.theta());
}
template<typename T> 
double disvars<T>::eprimeEnergy(){
  return peOutV4.e();
}
template<typename T> 
double disvars<T>::eprimepT(){
  return peOutV4.pT();
}
template<typename T> 
double disvars<T>::jetEnergy(){
  return jet4vec.e();
}
template<typename T> 
double disvars<T>::jetpT(){
  return jet4vec.pT();
}
template<typename T> 
double disvars<T>::jetTheta(){
  return 180/TMath::Pi()*(jet4vec.theta());
}
template<typename T> 
double disvars<T>::xBjorken(){
  return Q2()/(2.*pProtonV4*pPhoton);
}
template<typename T> 
double disvars<T>::yLine(){
  return (pProtonV4*pPhoton)/(pProtonV4*peInV4);
}
template<typename T> 
bool isBHadron(int id){
  //This snippet is menat to capture all B hadrons
  //as given in the PDG
  if(id<0) id*=-1;
  if(id<500) return false;
  return (fmod(id/100,5.)==0.0 || id/1000==5);
}
template<typename T> 
double disvars<T>::myRapidity(Vec4 p){
  return 0.5*log(p.pPos()/p.pNeg());
}
#endif
