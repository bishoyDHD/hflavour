#ifndef storage_h
#define storage_h 1
// This is merely a file for storing TTree kinematic vars

typedef struct{
  Double_t Q2,W2;
  Double_t yline,x_Bjor;
  Double_t E_e,theta_e,pT_e,rapidity_e;
  // jet kinematics
  std::vector<Double_t> theta_jet,E_jet;
  std::vector<Double_t> rapidity_jet,pT_jet,pTdist_jet;
  // parton kinematics
  std::vector<Double_t> theta_parton,E_parton;
  std::vector<Double_t> rapidity_parton,pT_parton,pTdist_parton;

  // Heavy quark kinematics
  std::vector<Double_t> pT_bquark,pT_bbarquark;
  std::vector<Double_t> pT_B0meson,pT_B0barmeson;
  std::vector<Double_t> pT_Bplusmeson,pT_Bminusmeson;
  std::vector<Double_t> pTCDFrap_Bplusmeson,pTCDFrap_Bminusmeson;
  std::vector<Double_t> pT_b2eminus,pT_b2eplus; // e-,e+ from b
} diskinematics;

#endif
