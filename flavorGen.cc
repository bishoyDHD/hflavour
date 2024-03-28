// flavorGen.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Bishoy (DHD)
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; electron-proton; event shapes; heavy meson production,

// This is a simple test program.
// DIS for EIC physics and study of Heavy Flavor generation via gamma-gluon fusion.

#include "Pythia8/Pythia.h"
#include "Pythia8/Basics.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TH1D.h>
#include "disvars.h"
#include "storage.h"

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;
  Event& event = pythia.event;

  double E_proton=250.;
  double E_electron=20.0;
  double Q2_min=0.;
  // set default no. of events
  int nEvents=1500000;

  //pythia.readString("Random.seed = 3000000");

  // Set up incoming beam for frame with unequal beam energies
  pythia.readString("Beams:frameType=2");
  // Proton beam config
  pythia.readString("Beams:idA = 2212");
  pythia.settings.parm("Beams:eA",E_proton);
  // e- beam config
  pythia.readString("Beams:idB =  11");
  pythia.settings.parm("Beams:eB",E_electron);

  // DIS Process selection.
  // Turn on only bbar production:
  // g g    -> b bbar (subprocess 123)
  // q qbar -> b bbar (subprocess 124)
  pythia.readString("HardQCD:all = on");
  //Charged current
  pythia.readString("WeakBosonExchange:ff2ff(t:W)=on");
  // Neutral current (with gamma/Z interference).
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  // phase-space cut: minimal Q2 of process
  pythia.settings.parm("PhaseSpace:Q2Min",Q2_min);
  pythia.readString("23:mMin=1.0");
 
  // Random number Generator Should be Set Here if needed (before pythia.init())
  // On seeds:
  // seed = -1 : default (any number < 0 will revert to the default).  seed = 19780503
  // seed = 0 : calls Stdlib time(0) to provide a seed based on the unix time
  // seed = 1 through 900 000 000: various numbers that can be used as seeds
 
  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");
  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");
  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // Initialize()
  pythia.init();

  // create TTree and branches containing DIS kinematics
  diskinematics diskin;
  TTree *ftree=new TTree("hflavor","EIC DIS Kinematics");
  //ftree->SetAutoSave();
  ftree->Branch("Q2",&diskin.Q2);
  ftree->Branch("W2",&diskin.W2);
  ftree->Branch("yline",&diskin.yline);
  ftree->Branch("x_Bjor",&diskin.x_Bjor);
  ftree->Branch("E_e",&diskin.E_e);
  ftree->Branch("theta_e",&diskin.theta_e);
  ftree->Branch("rapidity_e",&diskin.rapidity_e);
  ftree->Branch("pT_e",&diskin.pT_e);
  // jet kinematics
  ftree->Branch("E_jet",&diskin.E_jet);
  ftree->Branch("pT_jet",&diskin.pT_jet);
  ftree->Branch("pTdist_jet",&diskin.pTdist_jet);
  ftree->Branch("theta_jet",&diskin.theta_jet);
  ftree->Branch("rapidity_jet",&diskin.rapidity_jet);
  // parton kinematics
  ftree->Branch("E_parton",&diskin.E_parton);
  ftree->Branch("pT_parton",&diskin.pT_parton);
  ftree->Branch("pTdist_parton",&diskin.pTdist_parton);
  ftree->Branch("theta_parton",&diskin.theta_parton);
  ftree->Branch("rapidity_parton",&diskin.rapidity_parton);
  // heavy flavor kinematics
  ftree->Branch("pT_bquark",&diskin.pT_bquark);
  ftree->Branch("pT_bbarquark",&diskin.pT_bbarquark);
  // Histograms.
  Hist mult("charged multiplicity", 100, -0.5, 799.5);
  double W_max=std::sqrt(4.*E_proton*E_electron);
  // Set up ROOT file & Histos
  //===============================================================
  //      x-Q2 kinematics for scattered particles
  //===============================================================
  TH1D* h1Q = new TH1D("Q [GeV]","",100,-0.5,49.5);
  TH1D* h1W = new TH1D("W [GeV]","",100,0.,W_max);
  TH1D* h1x = new TH1D("h1x","",100,0.,1);
  TH1D* h1y = new TH1D("h1y","",100,0.,1);
  TH1D* h1pTe = new TH1D("h1pTe","",100,0.,20.);
  TH1D* h1pTparton = new TH1D("h1pTparton","",100,0.,50.);
  TH1D* h1pTdist = new TH1D("h1pTdist","ratio pT_parton/pT_electron",100,0.,5.);

  double theta_eprime[]={5,50,90,100,115,120,170};
  double E_eprime[]={5,10,15,21,23,30,50,120,170};
  double theta_jet[]={10,30,50,90,120,150,170};
  double E_jet[]={5,20,25,50,90,120,170};
  // vars for angle and energy cuts
  double thetaVal=0.0,angLow=0.0, angUp=0.0;
  double energyVal=0.0,eprimeVal=0.0,eprimeUp=0.0,eprimeLow=0.0;
  double energyLow=0.0, energyUp=0.0;
  std::vector<std::vector<double>> x_eprimeE(7),Q2_eprimeE(7),x_Beprime(7),Q2_eprime(7);
  std::vector<std::vector<double>> x_jetE(7),Q2_jetE(7),jetE(7),x_Bjet(7),Q2_jet(7);

  // Begin event loop. Generate event. Skip if error. List first few.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;

    // empty jet vectors                                // clear parton vectors
    diskin.theta_jet.clear();                           diskin.theta_parton.clear();
    diskin.E_jet.clear();                               diskin.E_parton.clear();
    diskin.rapidity_jet.clear();                        diskin.rapidity_parton.clear();
    diskin.pT_jet.clear();                              diskin.pT_parton.clear();
    diskin.pTdist_jet.clear();                          diskin.pTdist_parton.clear(); 

    if (iEvent < 1) {pythia.info.list(); pythia.event.list();}
    
    //obtain DIS kinematics
    // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
    // pProton = event[1].p();
    // peIn    = event[2].p();
    // peOut   = event[6].p();
    disvars<Vec4> disEvts(event[1].p(),event[2].p(),event[6].p());
    thetaVal=disEvts.eprimeTheta();
    eprimeVal=disEvts.eprimeEnergy();
    h1Q->Fill(disEvts.Q2());
    h1W->Fill(disEvts.W2());
    h1x->Fill(disEvts.xBjorken());
    h1y->Fill(disEvts.yLine());
    h1pTe->Fill(disEvts.eprimepT());
    // Assign eprime tree variables
    diskin.Q2=disEvts.Q2();
    diskin.W2=disEvts.W2();
    diskin.x_Bjor=disEvts.xBjorken();
    diskin.yline=disEvts.yLine();
    diskin.E_e=disEvts.eprimeEnergy();
    diskin.theta_e=disEvts.eprimeTheta();
    diskin.pT_e=disEvts.eprimepT();
    diskin.rapidity_e=disEvts.myRapidity(event[2].p());

    for(int i = 0; i < event.size(); ++i){
      // scattered non-lepton (parton) kinematic information
      if(event[i].statusAbs() == 23 && event[i].isParton()){
        disEvts.setJet4Vec(event[i].p());
        // Assign parton tree variables
        diskin.theta_parton.push_back(disEvts.jetTheta());
        diskin.E_parton.push_back(disEvts.jetEnergy());
        diskin.rapidity_parton.push_back(disEvts.myRapidity(event[i].p()));
        diskin.pT_parton.push_back(disEvts.jetpT());
        diskin.pTdist_parton.push_back(disEvts.jetpT()/disEvts.eprimepT());
      }
      // scattered jet kinematic information
      if(event[i].statusAbs() == 43 && event[i].isDiquark()){
        disEvts.setJet4Vec(event[i].p());
        thetaVal=disEvts.jetTheta();
        energyVal=disEvts.jetEnergy();
        h1pTparton->Fill(disEvts.jetpT());
        h1pTdist->Fill(disEvts.jetpT()/disEvts.eprimepT());
        // Assign parton tree variables
        diskin.theta_jet.push_back(disEvts.jetTheta());
        diskin.E_jet.push_back(disEvts.jetEnergy());
        diskin.rapidity_jet.push_back(disEvts.myRapidity(event[i].p()));
        diskin.pT_jet.push_back(disEvts.jetpT());
        diskin.pTdist_jet.push_back(disEvts.jetpT()/disEvts.eprimepT());
      }
    }// end of parton kinematics for-loop
  // End of event loop. Statistics. Histogram. Done.
  ftree->Fill();
  }// event loop
  pythia.stat();

  //Write Output ROOT hisotgram into ROOT file
  TFile* outFile = new TFile("pythiaOutputHistos1M.root","RECREATE");

  h1Q->Write();
  h1W->Write();
  h1x->Write();
  h1y->Write();
  h1pTe->Write();
  h1pTparton->Write();
  h1pTdist->Write();
  ftree->Write();

  outFile->Close();

  // Done.
  return 0;
}
