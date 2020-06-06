#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"


void MakeCMSDataHisto() {

  //gStyle->SetErrorX(0);

  TFile* f1 = new TFile("../yskimPA1st_OpSign_20177262037_unIdentified.root");
  TTree* hf = (TTree*)f1->Get("hf");

  TFile* f2 = new TFile("../yskimPA2nd_OpSign_20177262044_unIdentified.root");
  TTree* hf2 = (TTree*)f2->Get("hf");

  TFile* outfile = new TFile("pPbActivityHistos.root","RECREATE");

  int maxNtracks = 330;
  int minNtracks = 0;
  float maxSumET = 270;
  float minSumET = 0.0;

  TH1F* hNtracks = new TH1F("hNtracks","hNtracks",165,minNtracks,maxNtracks);
  TH1F* hNtracks2 = new TH1F("hNtracks2","hNtracks2",165,minNtracks,maxNtracks);
  TH1F* hSumET_HF = new TH1F("hSumET_HF","hSumET_HF",100,minSumET,maxSumET);
  TH1F* hSumET_HF2 = new TH1F("hSumET_HF2","hSumET_HF2",100,minSumET,maxSumET);
  
  TCanvas* c1 = new TCanvas("cNtracks","cNtracks",4,45,550,520);
  c1->SetLogy();
  c1->cd();
  hf2->Draw("Ntracks>>hNtracks2","","E");
  hf->Draw("Ntracks>>hNtracks","","E");
  hNtracks->Add(hNtracks2);

  TCanvas* c2 = new TCanvas("cSumET","cSumET",604,45,550,520);
  c2->SetLogy();
  c2->cd();
  hf2->Draw("SumET_HF>>hSumET_HF2","","E");
  hf->Draw("SumET_HF>>hSumET_HF","","E");
  hSumET_HF->Add(hSumET_HF2);

  outfile->cd();
  hNtracks->Write();
  hSumET_HF->Write();
  //outfile->Close();
}
