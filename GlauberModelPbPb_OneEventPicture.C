#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TH2.h"
#include "time.h"
#include "TNtuple.h"
#include "TFile.h"

void GlauberModelPbPb_OneEventPicture(){

  gStyle->SetOptStat(0);

  //set constant parameters
  const double pi = TMath::Pi();
  const double R = 6.62;
  const double a = 0.546;
  const double max = 2*R;
  const double sigma = 7.0; //70 mb (1 b = 10^-28 m^2 = 100 fm^2)
  const double dtest = sigma/pi;
  const double scale = 5;
  const int numNucleons = 208;

  //create probability distributions
  TF1 *prob = new TF1("prob","(x**2)/(1+TMath::Exp((x-[0])/[1]))",0,max);
  prob->SetParameters(R,a);
  TF1 *thprob = new TF1("thprob","sin(x)",0,pi);
  TF1 *bprob = new TF1("bprob","x",0,1.5*max);

  //initialize variables
  double x1[numNucleons], y1[numNucleons], x2[numNucleons], y2[numNucleons];
  double r, theta, phi, dx, dy, b;
  TRandom3 rnd3(0);
  int numpart, numcol;
  int iEvents = 0;

  //b = bprob->GetRandom();
  b = 4;

  //populate nucleus with nucleons
  for (int i = 0; i<numNucleons; i++) {
    r = prob->GetRandom();
    theta = thprob->GetRandom();
    phi = 2*pi*rnd3.Rndm();
    x1[i] = r*TMath::Sin(theta)*TMath::Cos(phi);
    y1[i] = r*TMath::Sin(theta)*TMath::Sin(phi);
    r = prob->GetRandom();
    theta = thprob->GetRandom();
    phi = 2*pi*rnd3.Rndm();
    x2[i] = r*TMath::Sin(theta)*TMath::Cos(phi)+b;
    y2[i] = r*TMath::Sin(theta)*TMath::Sin(phi);
  }

  //search for collisions and count participants
  numcol = 0;
  numpart = 0;
  int participants1[numNucleons] = {0};
  int participants2[numNucleons] = {0};
  for (int i = 0; i<numNucleons; i++) {
    for (int j = 0; j<numNucleons; j++) {
      dx = x1[i]-x2[j];
      dy = y1[i]-y2[j];
      if (dx*dx + dy*dy < dtest) {
        numcol++;
        participants1[i] = 1;
        participants2[j] = 1;
      }
    }
  }
  for (int i = 0; i<numNucleons; i++) {
    numpart += participants1[i]+participants2[i];
  }

  //Make a pretty picture
  double histxmin = -max+b/2;
  double histxmax = max+b/2;
  double histymin = -max;
  double histymax = max;
  TH2D* nuc1 = new TH2D("nuc1","nuc1",100,histxmin,histxmax,100,histymin,histymax);
  TH2D* nuc2 = new TH2D("nuc2","nuc2",100,histxmin,histxmax,100,histymin,histymax);
  TH2D* parts = new TH2D("parts","parts",100,histxmin,histxmax,100,histymin,histymax);
  for (int i = 0; i<numNucleons; i++) {
    if (participants1[i]) parts->Fill(x1[i],y1[i]);
    else nuc1->Fill(x1[i],y1[i]);
    if (participants2[i]) parts->Fill(x2[i],y2[i]);
    else nuc2->Fill(x2[i],y2[i]);
  }

  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  nuc1->GetXaxis()->SetTitle("x (fm)");
  nuc1->GetYaxis()->SetTitle("y (fm)");
  nuc1->SetTitle("Glauber Model Event");
  nuc1->SetMarkerColor(8);
  nuc1->SetMarkerStyle(4);
  nuc1->SetMarkerSize(1.5);
  nuc1->Draw();
  nuc2->SetMarkerColor(9);
  nuc2->SetMarkerStyle(4);
  nuc2->SetMarkerSize(1.5);
  nuc2->Draw("same");
  parts->SetMarkerColor(2);
  parts->SetMarkerStyle(4);
  parts->SetMarkerSize(1.5);
  parts->Draw("same");

  TEllipse *el1 = new TEllipse(0,0,7,7);
  el1->SetFillStyle(4000);
  el1->SetLineStyle(10);
  el1->Draw();
  TEllipse *el2 = new TEllipse(b,0,7,7);
  el2->SetFillStyle(4000);
  el2->SetLineStyle(10);
  el2->Draw();

  double xpos = 0.05*(histxmax-histxmin)+histxmin;
  double ypos = 0.95*(histymax-histymin)+histymin;
  double dy = 0.06*(histymax-histymin);
  TText *t1 = new TText(xpos,ypos,Form("A = %i",numNucleons));
  t1->SetTextFont(43);
  t1->SetTextSize(20);
  t1->Draw();
  TLatex *t2 = new TLatex(xpos,ypos-dy,Form("N_{part} = %i",numpart));
  t2->SetTextFont(43);
  t2->SetTextSize(20);
  t2->Draw();
  TLatex *t3 = new TLatex(xpos,ypos-2*dy,Form("N_{coll} = %i",numcol));
  t3->SetTextFont(43);
  t3->SetTextSize(20);
  t3->Draw();

  c1->SaveAs(Form("GMEventPictureb%.0f.png",b));

}
