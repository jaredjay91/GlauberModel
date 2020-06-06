#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"


using namespace std;

void Make_nice_histos() {

  //import histograms
  TString histFileName = "GlauberModelpPbData100000000.root";
  cout << histFileName << endl;
  TFile *histFile = new TFile(histFileName);
  TNtuple* ntupleResults = (TNtuple*)histFile->Get("ntuple;1");

  const int min = 0;
  const int max = 30;

  //generate canvas with all 3 hisograms
  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(3,1);
  TH1F * histob = new TH1F("histob", "Impact Parameter", 100,min,15);
  TH1F * histonumpart = new TH1F("histonumpart", "Number of Participants", 30,min,max);
  TH1F * histonumcol = new TH1F("histonumcol", "Number of Collisions", 30,min,max);

  cntuple->cd(1);
  histob->SetXTitle("Impact Parameter");
  histob->GetXaxis()->SetTitleSize(0.05);
  histob->GetYaxis()->SetLabelSize(0.05);
  histob->GetXaxis()->SetLabelSize(0.05);
  histob->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("b>>histob");

  cntuple->cd(2);
  histonumpart->SetXTitle("Number of Participants");
  histonumpart->GetXaxis()->SetTitleSize(0.05);
  histonumpart->GetYaxis()->SetLabelSize(0.05);
  histonumpart->GetXaxis()->SetLabelSize(0.05);
  histonumpart->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("numpart>>histonumpart");

  cntuple->cd(3);
  histonumcol->SetXTitle("Number of Collisions");
  histonumcol->GetXaxis()->SetTitleSize(0.05);
  histonumcol->GetYaxis()->SetLabelSize(0.05);
  histonumcol->GetXaxis()->SetLabelSize(0.05);
  histonumcol->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("numcol>>histonumcol");

  cntuple->SaveAs("bNumpartNumcol.png");

}
