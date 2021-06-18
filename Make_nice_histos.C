#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"


using namespace std;

void Make_nice_histos() {

  //import histograms
  TString histFileName = "GlauberModelData1000000.root";
  cout << histFileName << endl;
  TFile *histFile = new TFile(histFileName);
  TNtuple* ntupleResults = (TNtuple*)histFile->Get("ntuple;1");


  //generate canvas with all 3 hisograms
  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(3,1);
  TH1F * histob = new TH1F("histob", "Impact Parameter", 100,0,21);
  TH1F * histonumpart = new TH1F("histonumpart", "Number of Participants", 100,0,500);
  TH1F * histonumcol = new TH1F("histonumcol", "Number of Collisions", 100,0,3000);

  cntuple->cd(1);
  histob->SetTitle("");
  histob->SetXTitle("b");
  histob->GetXaxis()->SetTitleSize(0.04);
  histob->GetYaxis()->SetLabelSize(0.04);
  histob->GetXaxis()->SetLabelSize(0.04);
  //histob->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("b>>histob");

  cntuple->cd(2);
  histonumpart->SetTitle("");
  histonumpart->SetXTitle("N_{part}");
  histonumpart->GetXaxis()->SetTitleSize(0.04);
  histonumpart->GetYaxis()->SetLabelSize(0.04);
  histonumpart->GetXaxis()->SetLabelSize(0.04);
  //histonumpart->GetXaxis()->SetRangeUser(min,max);
  gPad->SetLogy();
  ntupleResults->Draw("numpart>>histonumpart");

  cntuple->cd(3);
  histonumcol->SetTitle("");
  histonumcol->SetXTitle("N_{coll}");
  histonumcol->GetXaxis()->SetTitleSize(0.04);
  histonumcol->GetYaxis()->SetLabelSize(0.04);
  histonumcol->GetXaxis()->SetLabelSize(0.04);
  //histonumcol->GetXaxis()->SetRangeUser(min,max);
  gPad->SetLogy();
  ntupleResults->Draw("numcol>>histonumcol");

  cntuple->SaveAs("bNumpartNumcol.png");

}
