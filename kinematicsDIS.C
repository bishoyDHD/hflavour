#include <vector>
 
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"

void DrawTag(TGraph *g, int i, TString t){
   double x = g->GetPointX(i);
   double y = g->GetPointY(i);
   auto tag = new TLatex(x,y,t.Data());
   tag->Draw();
}

void kinematicsDIS(){
  double Q2,W2,yline,x_Bjor,E_e,theta_e,pT_e,rapidity_e;
  std::vector<Float_t> *theta_parton=nullptr,*E_parton=nullptr,*rapidity_parton=nullptr;
  std::vector<Float_t> *pT_parton=nullptr,*pTdist_parton=nullptr;
  //Jets
  std::vector<Float_t> *theta_jet=nullptr,*E_jet=nullptr,*rapidity_jet=nullptr;
  std::vector<Float_t> *pT_jet=nullptr,*pTdist_jet=nullptr;

  TFile *f1=new TFile("pythiaOutputHistos1M.root","read");
  TTree *t1 = (TTree*)f1->Get("hflavor");

  t1->SetBranchAddress("Q2",&Q2);
  t1->SetBranchAddress("W2",&W2);
  t1->SetBranchAddress("yline",&yline);
  t1->SetBranchAddress("x_Bjor",&x_Bjor);
  t1->SetBranchAddress("E_e",&E_e);
  t1->SetBranchAddress("theta_e",&theta_e);
  t1->SetBranchAddress("rapidity_e",&rapidity_e);
  t1->SetBranchAddress("pT_e",&pT_e);
  t1->SetBranchAddress("E_parton",&E_parton);
  t1->SetBranchAddress("pT_parton",&pT_parton);
  t1->SetBranchAddress("pTdist_parton",&pTdist_parton);
  t1->SetBranchAddress("theta_parton",&theta_parton);
  t1->SetBranchAddress("rapidity_parton",&rapidity_parton);
  //Jets
  t1->SetBranchAddress("E_jet",&E_jet);
  t1->SetBranchAddress("pT_jet",&pT_jet);
  t1->SetBranchAddress("pTdist_jet",&pTdist_jet);
  t1->SetBranchAddress("theta_jet",&theta_jet);
  t1->SetBranchAddress("rapidity_jet",&rapidity_jet);
  // setting up histograms
  TH1D* h1Q = new TH1D("Q [GeV]","",100,-0.5,49.5);
  //TH1D* h1W = new TH1D("W [GeV]","",100,0.,W_max);
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


  double theta_eprime[]={50,80,110,140,150,160,170,180};
  double E_eprime[]={5,10,15,21,23,30,50,120,170};
  double theta_prt[]={10,30,50,90,120,150,170};
  double E_prt[]={5,10,100,150,200,250,300,500};
  double W_line[]={5,10,20,50,75,80,100,125,150,190,200,270};
  double consty[]={1,0.5,0.1,0.01,0.001,0.0001,1e-5,1.0e-6};
  // vars for angle and energy cuts
  double thetaVal=0.0,angLow=0.0, angUp=0.0, yval=0.0,Wval=0.0;
  double energyVal=0.0,eprimeVal=0.0,eprimeUp=0.0,eprimeLow=0.0;
  double energyLow=0.0, energyUp=0.0;
  double EprtVal=0.0, thetaPrt=0.0;

  std::vector<std::vector<double>> x_eprimeE(7),Q2_eprimeE(7),x_Beprime(7),Q2_eprime(7),W(7);
  std::vector<std::vector<double>> x_yline(7),Q2_yline(7),W_yline(7),y_(7);
  std::vector<std::vector<double>> x_W(7),Q2_W(7);
  std::vector<std::vector<double>> x_jetE(7),Q2_jetE(7),jetE(7),x_Bjet(7),Q2_jet(7);
  std::vector<std::vector<double>> x_prtE(7),Q2_prtE(7),prtE(7),x_Bprt(7),Q2_prt(7);
  Int_t nentries=t1->GetEntries();
  for(int i=0; i<nentries; i++){
    t1->GetEntry(i);

    //store information for scatter angle
    thetaVal=theta_e;
    eprimeVal=E_e;
    yval=yline;
    Wval=std::sqrt(W2);
    for(int iangle=0; iangle<6; iangle++){
      angLow=theta_eprime[iangle]-0.5;
      angUp=theta_eprime[iangle]+0.5;
      if(thetaVal>=angLow && thetaVal<=angUp){
        x_Beprime[iangle].push_back(x_Bjor);
        Q2_eprime[iangle].push_back(Q2);
      }
      // Energy of scattered electron
      eprimeLow=E_eprime[iangle]-0.5;
      eprimeUp=E_eprime[iangle]+0.5;
      if(eprimeVal>=15 && eprimeVal<=22){
        eprimeLow=E_eprime[iangle]-0.05;
        eprimeUp=E_eprime[iangle]+0.05;
      }
      if(eprimeVal>=eprimeLow && eprimeVal<=eprimeUp){
        x_eprimeE[iangle].push_back(x_Bjor);
        Q2_eprimeE[iangle].push_back(Q2);
      }
      // lines of constant y
      if(yval>=consty[iangle]-0.0001 && yval<=consty[iangle]+0.0001){
        x_yline[iangle].push_back(x_Bjor);
        Q2_yline[iangle].push_back(Q2);
        W_yline[iangle].push_back(std::sqrt(W2));
        y_[iangle].push_back(yline);
      }
      // plots of W
      if(Wval>=W_line[iangle]-0.5 && Wval<=W_line[iangle]+0.5){
        x_W[iangle].push_back(x_Bjor);
        Q2_W[iangle].push_back(Q2);
      }
      // parton kinematics
      // Energy of scattered parton
      eprimeLow=E_prt[iangle]-0.005;
      eprimeUp=E_prt[iangle]+0.005;
      angLow=theta_prt[iangle]-0.005;
      angUp=theta_prt[iangle]+0.005;
      for(int nval=0; nval<E_parton->size(); nval++){
        EprtVal=E_parton->at(nval);
        thetaPrt=theta_parton->at(nval);
        //std::cout<<"Value at E_parton["<<nval<<"] = "<<EprtVal<<std::endl;
        if(EprtVal>=eprimeLow && EprtVal<=eprimeUp){
          x_prtE[iangle].push_back(x_Bjor);
          Q2_prtE[iangle].push_back(Q2);
        }
        if(thetaPrt>=angLow && thetaPrt<=angUp){
          x_Bprt[iangle].push_back(x_Bjor);
          Q2_prt[iangle].push_back(Q2);
        }
      }
      // Energy of scattered jet
      for(int nval=0; nval<E_jet->size(); nval++){
        EprtVal=E_jet->at(nval);
        //std::cout<<"Value at E_parton["<<nval<<"] = "<<EprtVal<<std::endl;
        if(EprtVal>=eprimeLow && EprtVal<=eprimeUp){
          x_jetE[iangle].push_back(x_Bjor);
          Q2_jetE[iangle].push_back(Q2);
        }
      }
    }// end of electron kinematics for-loop
  }

/*
  TCanvas* e_kin=new TCanvas("e_kin"," ",200,10,700,500);
  e_kin->cd();
  h1Q->Draw();
  e_kin->SetLogx();
  e_kin->SetLogy();
  e_kin->Update();
*/
  TString tagtheta_eprime[]={"30","50","100","150","160","170","200"};
/*
  TString tagE_eprime[]={5,10,15,21,23,30,50,120,170};
  TString tagtheta_prt[]={10,30,50,90,120,150,170};
  TString tagE_prt[]={5,10,100,150,200,250,300,500};
  TString tagW_line[]={5,10,20,50,75,80,100,125,150,190,200,270};
  TString tagconsty[]={1,0.5,0.1,0.01,0.001,0.0001,1e-5,1.0e-6};
*/
  TGraph* grTheta_eprime[10],*grE_eprime[10],*grTheta_jet[10],*grE_jet[10];
  TGraph* grTheta_prt[10],*grE_prt[10];
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
    n=Q2_prtE[i].size();
    grE_prt[i]=new TGraph(n,&(x_prtE[i][0]),&(Q2_prtE[i][0]));
    grE_prt[i]->SetMarkerStyle(kFullDotSmall);
    // ----- Parton Angles
    n=Q2_prt[i].size();
    grTheta_prt[i]=new TGraph(n,&(x_Bprt[i][0]),&(Q2_prt[i][0]));
    grTheta_prt[i]->SetMarkerStyle(kFullDotSmall);
    //Jet Graphs
    n=Q2_jetE[i].size();
    grE_jet[i]=new TGraph(n,&(x_jetE[i][0]),&(Q2_jetE[i][0]));
    grE_jet[i]->SetMarkerStyle(kFullDotSmall);
    //W-const yline
    n=Q2_W[i].size();
    grW[i]=new TGraph(n,&(x_W[i][0]),&(Q2_W[i][0]));
    grW[i]->SetMarkerStyle(kFullDotSmall);
  }

  TCanvas* ang = new TCanvas("ang","scattered electron Angle",700,900);
  ang->cd();
  TMultiGraph* mgAng = new TMultiGraph();
  auto legTheta_e = new TLegend(0.6,0.2,0.8,0.7);
  mgAng->SetTitle("Scattered Electron Angle");
  for(int i=0; i<6; i++){
    mgAng->Add(grTheta_eprime[i]);
    legTheta_e->AddEntry(grTheta_eprime[i],tagtheta_eprime[i],"LP");
    DrawTag(grTheta_eprime[i],i,tagtheta_eprime[i]);
  }
  ang->SetLogx();
  ang->SetLogy();
  mgAng->GetXaxis()->SetTitle("x");;
  mgAng->GetYaxis()->SetTitle("Q^{2} [GeV]^{2}");;
  mgAng->Draw("AP");
  //legTheta_e->Draw();
  //Parton -------------
  TCanvas* prt = new TCanvas("parton","Parton Energy",700,900);
  prt->cd();
  TMultiGraph* eprt = new TMultiGraph();
  eprt->SetTitle("Scattered Parton Energy");
  for(int i=0; i<6; i++)
    eprt->Add(grE_prt[i]);
  prt->SetLogx();
  prt->SetLogy();
  eprt->GetXaxis()->SetTitle("x");;
  eprt->GetYaxis()->SetTitle("Q^{2} [GeV]^{2}");;
  eprt->Draw("AP");

  TCanvas* prt_ang = new TCanvas("prtAng","Parton Angle",700,900);
  prt_ang->cd();
  TMultiGraph* angprt = new TMultiGraph();
  angprt->SetTitle("Scattered Parton Angle");
  for(int i=0; i<6; i++)
    angprt->Add(grTheta_prt[i]);
  prt_ang->SetLogx();
  prt_ang->SetLogy();
  angprt->GetXaxis()->SetTitle("x");;
  angprt->GetYaxis()->SetTitle("Q^{2} [GeV]^{2}");;
  angprt->Draw("AP");
  //Jet
  TCanvas* jet = new TCanvas("jet","Jet Energy",700,900);
  jet->cd();
  TMultiGraph* ejet = new TMultiGraph();
  ejet->SetTitle("Jet Energy");
  for(int i=0; i<4; i++)
    ejet->Add(grE_jet[i]);
  jet->SetLogx();
  jet->SetLogy();
  ejet->Draw("AP");

  TCanvas* ene = new TCanvas("ene","Scattered electron Energy",700,900);
  ene->cd();
  TMultiGraph* mgEne = new TMultiGraph();
  mgEne->SetTitle("Scattered Electron Energy");
  for(int i=0; i<6; i++)
    mgEne->Add(grE_eprime[i]);
  ene->SetLogx();
  ene->SetLogy();
  mgEne->GetXaxis()->SetTitle("x");;
  mgEne->GetYaxis()->SetTitle("Q^{2} [GeV]^{2}");;
  mgEne->Draw("AP");

  TCanvas* wy = new TCanvas("wy","W and Constant y",700,900);
  wy->cd();
  TMultiGraph* mgWy = new TMultiGraph();
  mgWy->SetTitle("Lines of Constant #font[12]{y}");
  for(int i=0; i<6; i++)
    mgWy->Add(grW[i]);
  wy->SetLogx();
  wy->SetLogy();
  mgWy->GetXaxis()->SetTitle("x");;
  mgWy->GetYaxis()->SetTitle("Q^{2} [GeV]^{2}");;
  mgWy->Draw("AP");

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

}
