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

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;
  Event& event = pythia.event;

  double E_proton=250.;
  double E_electron=18.0;
  double Q2_min=0.;
  // set default no. of events
  int nEvents=5000000;

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

  double Q2,W2,yline,x_Bjor,E_e,theta_e,pT_e,rapidity_e;
  std::vector<Float_t> theta_parton,E_jets,rapidity_jets,pT_parton,pTdist,etest;

  TTree *ftree=new TTree("hflavor","EIC DIS Kinematics");
  //ftree->SetAutoSave();
  ftree->Branch("Q2",&Q2);
  ftree->Branch("etest",&etest);
  ftree->Branch("W2",&W2);
  ftree->Branch("yline",&yline);
  ftree->Branch("x_Bjor",&x_Bjor);
  ftree->Branch("E_e",&E_e);
  ftree->Branch("theta_e",&theta_e);
  ftree->Branch("rapidity_e",&rapidity_e);
  ftree->Branch("pT_e",&pT_e);
  ftree->Branch("E_jets",&E_jets);
  ftree->Branch("pT_parton",&pT_parton);
  ftree->Branch("pTdist",&pTdist);
  ftree->Branch("theta_parton",&theta_parton);
  ftree->Branch("rapidity_jets",&rapidity_jets);
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

  /* Uncomment lines below for heavy flavor studies
  //===============================================================
  //      Heavy quark production from gamma-gluon fusion
  //===============================================================
  TH1D* multHist = new TH1D("multHist","Multiplicity",100,-0.5,99.5);
  TH1D* bquarkPt = new TH1D("bquarkPt","bquarkPt",100,0,50);
  TH1D* bbarquarkPt = new TH1D("bbarquarkPt","bbar quark Pt",100,0,50);
  TH1D* B0mesonPt = new TH1D("BOmesonPt","B0mesonPt",100,0,50);
  TH1D* B0barmesonPt = new TH1D("BObarmesonPt","B0bar meson Pt",100,0,50);
  TH1D* BplusmesonPt = new TH1D("BplusmesonPt","BplusmesonPt",100,0,50);
  TH1D* BminusmesonPt = new TH1D("BminusmesonPt","Bminus meson Pt",100,0,50);
  TH1D* BplusmesonPtCDFrap = new TH1D("BplusmesonPtCDFrap","BplusmesonPt |y|<1",100,0,50);
  TH1D* BminusmesonPtCDFrap = new TH1D("BminusmesonPtCDFrap","Bminus meson Pt |y|<1",100,0,50);
  TH1D* electronFrombPt = new TH1D("electronFrombPt","electrons from b",100,0,30);
  TH1D* positronFrombPt = new TH1D("positronFrombPt","positrons from b",100,0,30);
  TH1D* epluseminusMinv = new TH1D("epluseminusMinv","e+ e- Inv. Mass",100,0,30);
  TH1D* epluseminusRapidity = new TH1D("epluseminusRapidity","e+ e- y",80,-4,4);
  TH1D* epluseminusMinvMidRap = new TH1D("epluseminusMinvMidRap","e+ e- Inv. Mass |y|<0.5",300,0,30);
*/

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

    // empty vectors
    theta_parton.clear();
    E_jets.clear();
    rapidity_jets.clear();
    pT_parton.clear();
    pTdist.clear(); 
    etest.clear();

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
    Q2=disEvts.Q2();
    W2=disEvts.W2();
    x_Bjor=disEvts.xBjorken();
    yline=disEvts.yLine();
    E_e=disEvts.eprimeEnergy();
    theta_e=disEvts.eprimeTheta();
    pT_e=disEvts.eprimepT();
    rapidity_e=disEvts.myRapidity(event[2].p());
    etest.push_back(disEvts.eprimeEnergy());

    //store information for scatter angle
    for(int iangle=0; iangle<6; iangle++){
      angLow=theta_eprime[iangle]-0.5;
      angUp=theta_eprime[iangle]+0.5;
      if(thetaVal>=angLow && thetaVal<=angUp){
        x_Beprime[iangle].push_back(disEvts.xBjorken());
        Q2_eprime[iangle].push_back(disEvts.Q2());
      }
      // Energy of scattered electron
      eprimeLow=E_eprime[iangle]-0.5;
      eprimeUp=E_eprime[iangle]+0.5;
      if(eprimeVal>=15 && eprimeVal<=22){
        eprimeLow=E_eprime[iangle]-0.05;
        eprimeUp=E_eprime[iangle]+0.05;
      }
      if(eprimeVal>=eprimeLow && eprimeVal<=eprimeUp){
        x_eprimeE[iangle].push_back(disEvts.xBjorken());
        Q2_eprimeE[iangle].push_back(disEvts.Q2());
      }
    }// end of electron kinematics for-loop
    // scattered Jet kinematic information
    for(int i = 0; i < event.size(); ++i){
      if(event[i].statusAbs() == 43){
        disEvts.setJet4Vec(event[i].p());
        thetaVal=disEvts.jetTheta();
        energyVal=disEvts.jetEnergy();
        h1pTparton->Fill(disEvts.jetpT());
        h1pTdist->Fill(disEvts.jetpT()/disEvts.eprimepT());
        // Assign parton tree variables
        theta_parton.push_back(disEvts.jetTheta());
        E_jets.push_back(disEvts.jetEnergy());
        rapidity_jets.push_back(disEvts.myRapidity(event[i].p()));
        pT_parton.push_back(disEvts.jetpT());
        pTdist.push_back(disEvts.jetpT()/disEvts.eprimepT());
        for(int iangle=0; iangle<4; iangle++){
          angLow=theta_jet[iangle]-0.005;
          angUp=theta_jet[iangle]+0.005;
          if(thetaVal>=angLow && thetaVal<=angUp){
            x_Bjet[iangle].push_back(disEvts.xBjorken());
            Q2_jet[iangle].push_back(disEvts.Q2());
          }
          // Energy of current jet
          energyLow=E_jet[iangle]-0.5;
          energyUp=E_jet[iangle]+0.5;
          if(energyVal>=energyLow && energyVal<=energyUp){
            x_jetE[iangle].push_back(disEvts.xBjorken());
            Q2_jetE[iangle].push_back(disEvts.Q2());
          }
        }// end of parton-plotting for-loop
      }
    }// end of parton kinematics for-loop
  // End of event loop. Statistics. Histogram. Done.
  ftree->Fill();
  }// event loop
  pythia.stat();

  //Write Output ROOT hisotgram into ROOT file
  TFile* outFile = new TFile("pythiaOutputHistos1M.root","RECREATE");
  TGraph* grTheta_eprime[10],*grE_eprime[10],*grTheta_jet[10],*grE_jet[10];
  TGraph* grY[10],*grW[10];
  int n=0;
  for(int i=0; i<6; i++){
    n=Q2_eprime[i].size();
    grTheta_eprime[i]=new TGraph(n,&(x_Beprime[i][0]),&(Q2_eprime[i][0]));
    grTheta_eprime[i]->SetMarkerStyle(kFullDotSmall);
    // Energy
    n=Q2_eprimeE[i].size();
    grE_eprime[i]=new TGraph(n,&(x_eprimeE[i][0]),&(Q2_eprimeE[i][0]));
    grE_eprime[i]->SetMarkerStyle(kFullDotSmall);
    //Parton Graphs
    n=Q2_jet[i].size();
    grTheta_jet[i]=new TGraph(n,&(x_Bjet[i][0]),&(Q2_jet[i][0]));
    grTheta_jet[i]->SetMarkerStyle(kFullDotSmall);
  }

  TCanvas* ang = new TCanvas("ang","scattered electron Angle",700,900);
  ang->cd();
  TMultiGraph* mgAng = new TMultiGraph();
  for(int i=0; i<6; i++)
    mgAng->Add(grTheta_eprime[i]);
  ang->SetLogx();
  ang->SetLogy();
  mgAng->GetXaxis()->SetTitle("x");;
  mgAng->GetYaxis()->SetTitle("Q^{2}");;
  mgAng->Draw("AP");
  ang->Write();
  //Parton
  TCanvas* ang2 = new TCanvas("parton","",700,900);
  ang2->cd();
  TMultiGraph* mgAng2 = new TMultiGraph();
  for(int i=0; i<4; i++)
    mgAng2->Add(grTheta_jet[i]);
  ang2->SetLogx();
  ang2->SetLogy();
  mgAng2->Draw("AP");
  ang2->Write();

  TCanvas* ene = new TCanvas("ene","Scattered electron Energy",700,900);
  ene->cd();
  TMultiGraph* mgEne = new TMultiGraph();
  for(int i=0; i<6; i++)
    mgEne->Add(grE_eprime[i]);
  ene->SetLogx();
  ene->SetLogy();
  mgEne->GetXaxis()->SetTitle("x");;
  mgEne->GetYaxis()->SetTitle("Q^{2}");;
  mgEne->Draw("AP");
  ene->Write();
/*
  TCanvas* eprime = new TCanvas("eprime","",700,900);
  eprime->cd();
  TMultiGraph* epene = new TMultiGraph();
  epene->Add(ge_pe05);
  epene->Add(ge_pe50);
  epene->Add(ge_pe120);
  epene->Draw("AP");
  eprime->SetLogx();
  eprime->SetLogy();
  eprime->Write();*/
  h1Q->Write();
  h1W->Write();
  h1x->Write();
  h1y->Write();
  h1pTe->Write();
  h1pTparton->Write();
  h1pTdist->Write();
  ftree->Write();
/*
  h1pTparton->Write();
  multHist->Write();
  bquarkPt->Write();
  bbarquarkPt->Write();
  B0mesonPt->Write();
  B0barmesonPt->Write();
  BminusmesonPt->Write();
  BplusmesonPtCDFrap->Write();
  BminusmesonPtCDFrap->Write();
  electronFrombPt->Write();
  positronFrombPt->Write();
  epluseminusMinv->Write();
  epluseminusRapidity->Write();
  epluseminusMinvMidRap->Write();
*/
  outFile->Close();

  // Done.
  return 0;
}
