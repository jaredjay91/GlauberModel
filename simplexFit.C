#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "time.h"
#include "TFile.h"
#include "TText.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TLeaf.h"
#include "TLegend.h"

void simplexFit() {

  cout << "The program started." << endl;

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  gStyle->SetOptStat(0);

//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 1: GENERATE EVENTS WITH GLAUBER MODEL MONTE CARLO


  cout << "Loading Glauber Model Data." << endl;

  TFile *GMdata = new TFile("GlauberModelDataPbPb502TeV1000000.root");
  TNtuple* GMntuple = (TNtuple*)GMdata->Get("ntuple;1");
  TLeaf *bLeaf = GMntuple->GetLeaf("b");
  TLeaf *numpartLeaf = GMntuple->GetLeaf("numpart");
  TLeaf *numcolLeaf = GMntuple->GetLeaf("numcol");
  const int NEvents = (int)GMntuple->GetEntries();

  cout << "NEvents = " << NEvents << endl;

  cout << "Glauber Model loaded." << endl;


//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 2: APPLY PARTICLE PRODUCTION MODEL AND FIT MONTE CARLO TO CMS DATA

  //Create canvas
  TCanvas * c1 = new TCanvas("c1","c1",500,0,500,400);
  c1->SetLogy();
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.2);

  //Import CMS data.
  const char *inputFile = "10-001_Sum_ET_HF_Energy.root";
  //const char *inputFile = "HF2_PbPb_502TeV_MinBias_fromCrab_2020_11_03.root";
  TFile *thefile = new TFile(inputFile);
  TH1D* hNtracks_CMSdata = (TH1D*)thefile->Get("hSumEnergy_nominal");
  //TH1D* hNtracks_CMSdata = (TH1D*)thefile->Get("HF2Hist");
  double Ndata = hNtracks_CMSdata->GetEntries();
  int minBinFit = 5;
  int nHFbins = hNtracks_CMSdata->GetNbinsX();
  cout << "Number of bins = " << nHFbins << endl;

  //set poisson error bars in each bin
  hNtracks_CMSdata->Scale(Ndata);
  for (int i=1; i<nHFbins; i++){
    double poisErr = sqrt(hNtracks_CMSdata->GetBinContent(i));
    hNtracks_CMSdata->SetBinError(i,poisErr);
  }
  hNtracks_CMSdata->Sumw2();
  //hNtracks_CMSdata->Scale(1.0/Ndata);
  hNtracks_CMSdata->Scale(1.0/hNtracks_CMSdata->Integral());
  hNtracks_CMSdata->Draw("E");
  hNtracks_CMSdata->GetYaxis()->SetTitle("Count");
  hNtracks_CMSdata->GetXaxis()->SetTitle("SumET_HF");
  //hNtracks_CMSdata->SetMaximum(1e7);
  //hNtracks_CMSdata->SetMinimum(0.3);
  //hNtracks_CMSdata->SetMaximum(2e6);
  //hNtracks_CMSdata->SetMinimum(1e5);
  c1->SetBottomMargin(0.15);
  //c1->SetLeftMargin(0.2);

  cout << "total integral = " << hNtracks_CMSdata->Integral() << endl;

  cout << "Applying particle production model." << endl;

  //Define negative binomial distribution (NBD) with initial parameters
  //double mu = 1.43083;
  //double k = 0.753833;
  //double normval = 1.23;
  double mu = 1.5;
  double k = 1.0;
  double normval = 1.1975;
  double startingWidth = 0.2;
  double xscale = 1.0/1000;
  //double xscale = 1.0;
  //double xshift = 0;
  TF1 *NBD = new TF1("NBD","TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])) * (([0]/[1])**x)/(([0]/[1]+1)**(x+[1]))",0,28);
  NBD->SetParameters(mu,k);

  //Apply particle production model to get the SumEnergy histogram
  double SumEnergy;
  int numcol;
  double histomax = 5;
  //TH1D* hSumEnergy_100bins = new TH1D ("HSumEnergy100bins","Energy distribution in PbPb GM",165,0,histomax);
  //hSumEnergy_100bins->GetYaxis()->SetTitle("Count");
  TH1D* hSumEnergy_100bins = (TH1D*)hNtracks_CMSdata->Clone("HSumEnergy100bins");
  hSumEnergy_100bins->SetLineColor(2);


//-----------------------------------------------------
//DO THE FIT

  //Define arrays
  double Expected, Observed;
  double kGuess = k;
  double muGuess = mu;
  double kParam[3] = {kGuess*startingWidth, kGuess, kGuess*(1+startingWidth)}; //3 initial points
  double muParam[3] = {muGuess*startingWidth, muGuess*(1+startingWidth), muGuess};
  double normParam[3] = {normval, normval, normval};
  double chi2min = DBL_MAX;
  double chi2Val[3] = {chi2min, chi2min, chi2min};
  int smallOne, bigOne, otherOne;
  double dev = (kParam[2]-kParam[0])*(kParam[2]-kParam[0]) + (muParam[2]-muParam[0])*(muParam[2]-muParam[0]);

  //Use converging triangle method to find minimum negative log likelihood
  TFile* outFile = new TFile("fittedVals.root","recreate");
  TNtuple* ntuple = new TNtuple("ntuple","list of fitted values","k:mu:norm:chi2",1);
  TNtuple* ntuple1 = new TNtuple("ntuple1","list of fitted values","k:mu:norm:chi2",1);
  TNtuple* ntuple2 = new TNtuple("ntuple2","list of fitted values","k:mu:norm:chi2",1);
  TNtuple* ntuple3 = new TNtuple("ntuple3","list of fitted values","k:mu:norm:chi2",1);
  TNtuple* ntupleworst = new TNtuple("ntupleworst","list of fitted values","k:mu:norm:chi2",1);

  double pinch = 0.9;
  double maxOf3 = 0;
  double minOf3 = 0;

  int normN = 40;
  double normmin = 0.8;
  double normmax = 1.4;
  double normeps = (normmax-normmin)/normN;
  if (normN==1) {
    normmin = normval;
    normeps = 0.0;
  }

  while (dev>0.00001) {
  //while (dev>0.0001) {
    int changeOne = bigOne;
    smallOne = 0;
    otherOne = 0;
    bigOne = 0;
    maxOf3 = 0;
    minOf3 = chi2min;
    NBD->SetParameters(muParam[changeOne],kParam[changeOne]);
    //delete histogram and recreate it.
    delete hSumEnergy_100bins;
    hSumEnergy_100bins = (TH1D*)hNtracks_CMSdata->Clone("HSumEnergy100bins");
    for (int ipp = 0; ipp<NEvents; ipp++) {
      //For each nucleon collison (Ncol), sample from the NBD to estimate SumEnergy
      SumEnergy = 0;
      GMntuple->GetEntry(ipp);
      numcol = (int)numcolLeaf->GetValue();
      for (int icolcount = 0; icolcount<numcol; icolcount++) {
        SumEnergy += NBD->GetRandom()*xscale;
      }
      //cout << "SumEnergy = " << SumEnergy << endl;
      hSumEnergy_100bins->Fill(SumEnergy);
    }

    double chi2 = 0;
    double norm = normmin-normeps/2;
    double chi2minNorm = DBL_MAX;

    for (int inorm = 0; inorm < normN; inorm++) {
      norm += normeps;
      TH1D* hSumEnergy_norm = (TH1D*)hSumEnergy_100bins->Clone();
      hSumEnergy_norm->Scale(norm/NEvents);
      //cout << "total integral = " << hSumEnergy_100bins->Integral() << endl;
      chi2 = 0;

      //Calculate chi^2 for current choice of k, <n>.
      double Observed, Expected, diff, error;
      for (int ibin = minBinFit; ibin<nHFbins; ibin++) {
        Observed = (double)hNtracks_CMSdata->GetBinContent(ibin);
        Expected = (double)hSumEnergy_norm->GetBinContent(ibin);
        error = (double)hNtracks_CMSdata->GetBinError(ibin);
        if (error > 0) {
          diff = Observed - Expected;
          chi2 += diff*diff/(error*error);
          //cout << "Expected = " << Expected << endl;
        }
      }
      delete hSumEnergy_norm;
      if (chi2 < chi2minNorm) {
        chi2minNorm = chi2;
        normParam[changeOne] = norm;
      }
      chi2Val[changeOne] = chi2minNorm;
    }

    cout << "k = " << kParam[changeOne] << "; mu = " << muParam[changeOne] << "; norm = " << normParam[changeOne] << "; Chi^2 = " << chi2Val[changeOne] << endl;

    for (int j = 0; j<3; j++) {
      if (chi2Val[j]<=minOf3) {
        smallOne = j;
        minOf3 = chi2Val[j];
      }
      if (chi2Val[j]>=maxOf3) {
        bigOne = j;
        maxOf3 = chi2Val[j];
      }
      //cout << "(b,s) = " << "(" << muParam[j] << "," << kParam[j] << "): ";
      //cout << "chi2 = " << chi2[j] << endl;
    }

    otherOne = 3-smallOne-bigOne;

    if (minOf3<chi2min) chi2min = minOf3;
    //cout << Form("%i,%i,%i",smallOne,otherOne,bigOne) << endl;
    double x1 = muParam[smallOne];
    double y1 = kParam[smallOne];
    double x2 = muParam[otherOne];
    double y2 = kParam[otherOne];
    double x3 = muParam[bigOne];
    double y3 = kParam[bigOne];
    double xnew = x3+pinch*(x1+x2-2*x3);
    double ynew = y3+pinch*(y1+y2-2*y3);
    dev = (y3-y1)*(y3-y1) + (x3-x1)*(x3-x1);
    cout << "bigOne = " << bigOne << endl;
    cout << "smallOne = " << smallOne << endl;
    cout << "minOf3 = " << minOf3 << endl;
    cout << "maxOf3 = " << maxOf3 << endl;
    cout << "dev = " << dev << endl;
    muParam[bigOne] = xnew;
    kParam[bigOne] = ynew;
    //cout << "s = " << kParam[smallOne] << endl;
    //cout << "b = " << muParam[smallOne] << endl;
    //cout << "chi2min = " << chi2min << endl;
    ntuple->Fill(kParam[smallOne],muParam[smallOne],normParam[smallOne],minOf3);
    ntuple1->Fill(kParam[0],muParam[0],normParam[0],chi2Val[0]);
    ntuple2->Fill(kParam[1],muParam[1],normParam[1],chi2Val[1]);
    ntuple3->Fill(kParam[2],muParam[2],normParam[2],chi2Val[2]);
    ntupleworst->Fill(kParam[bigOne],muParam[bigOne],normParam[bigOne],maxOf3);
  }//end of while loop

  double kFitted = kParam[smallOne];
  double muFitted = muParam[smallOne];
  double normFitted = normParam[smallOne];
  cout << "(k,mu,norm) = ( " << kFitted << ", " << muFitted << ", " << normFitted << " )" << endl;
  cout << "chi2 = " << minOf3 << endl;


  //Draw path of parameters
  int Nsteps = (int)ntuple->GetEntries();
  const int constNsteps = Nsteps;
  float svals[constNsteps], bvals[constNsteps];
  float svals1[constNsteps], bvals1[constNsteps];
  float svals2[constNsteps], bvals2[constNsteps];
  float svals3[constNsteps], bvals3[constNsteps];
  float chi2bestvals[constNsteps];
  float chi2worstvals[constNsteps];
  float stepvals[constNsteps];
  cout << Nsteps << endl;
  TLeaf *kLeaf = ntuple->GetLeaf("k");
  TLeaf *muLeaf = ntuple->GetLeaf("mu");
  TLeaf *chi2bestLeaf = ntuple->GetLeaf("chi2");
  TLeaf *kLeaf1 = ntuple1->GetLeaf("k");
  TLeaf *muLeaf1 = ntuple1->GetLeaf("mu");
  TLeaf *kLeaf2 = ntuple2->GetLeaf("k");
  TLeaf *muLeaf2 = ntuple2->GetLeaf("mu");
  TLeaf *kLeaf3 = ntuple3->GetLeaf("k");
  TLeaf *muLeaf3 = ntuple3->GetLeaf("mu");
  TLeaf *chi2worstLeaf = ntupleworst->GetLeaf("chi2");
  for (int istep=0; istep<Nsteps; istep++) {
    ntuple->GetEntry(istep);
    svals[istep] = (float)kLeaf->GetValue();
    bvals[istep] = (float)muLeaf->GetValue();
    chi2bestvals[istep] = (float)chi2bestLeaf->GetValue();
    ntuple1->GetEntry(istep);
    svals1[istep] = (float)kLeaf1->GetValue();
    bvals1[istep] = (float)muLeaf1->GetValue();
    ntuple2->GetEntry(istep);
    svals2[istep] = (float)kLeaf2->GetValue();
    bvals2[istep] = (float)muLeaf2->GetValue();
    ntuple3->GetEntry(istep);
    svals3[istep] = (float)kLeaf3->GetValue();
    bvals3[istep] = (float)muLeaf3->GetValue();
    ntupleworst->GetEntry(istep);
    chi2worstvals[istep] = (float)chi2worstLeaf->GetValue();
    if (isinf(chi2worstvals[istep])) {
      chi2worstvals[istep] = 1e5;
    }
    cout << "chi2worstvals[" << istep << "] = " << chi2worstvals[istep] << endl;
    stepvals[istep] = (float)istep;
    cout << "stepvals[" << istep << "] = " << stepvals[istep] << endl;
  }

  TCanvas * c2 = new TCanvas("c2","c2",0,0,400,400);
  c2->cd();
  TGraph *gmin = new TGraph(Nsteps,bvals,svals);
  c2->SetLeftMargin(0.12);
  gmin->SetTitle("Path to minimum #chi^2");
  gmin->GetXaxis()->SetTitleSize(0.04);
  gmin->GetXaxis()->SetLabelSize(0.04);
  gmin->GetYaxis()->SetTitleSize(0.04);
  gmin->GetYaxis()->SetLabelSize(0.04);
  gmin->GetXaxis()->SetTitle("Parameter #mu");
  //gmin->GetXaxis()->SetTitleOffset(1.2);
  gmin->GetYaxis()->SetTitle("Parameter k");
  gmin->GetYaxis()->SetTitleOffset(1.5);
  //gmin->GetXaxis()->SetLimits(0.0,2*muFitted);
  //gmin->GetHistogram()->SetMinimum(0.0);
  //gmin->GetHistogram()->SetMaximum(2*kFitted);
  gmin->Draw("AC*");

  TCanvas * c3 = new TCanvas("c3","c3",0,400,400,400);
  c3->cd();
  c3->SetLeftMargin(0.12);
  TGraph *g1 = new TGraph(Nsteps,bvals1,svals1);
  g1->SetLineColor(kRed);
  TGraph *g2 = new TGraph(Nsteps,bvals2,svals2);
  g2->SetLineColor(kBlue);
  TGraph *g3 = new TGraph(Nsteps,bvals3,svals3);
  g3->SetLineColor(kGreen);
  g1->SetTitle("Convergence of the 3 points");
  g1->GetXaxis()->SetTitleSize(0.04);
  g1->GetXaxis()->SetLabelSize(0.04);
  g1->GetYaxis()->SetTitleSize(0.04);
  g1->GetYaxis()->SetLabelSize(0.04);
  g1->GetXaxis()->SetTitle("Parameter #mu");
  //g1->GetXaxis()->SetTitleOffset(1.2);
  g1->GetYaxis()->SetTitle("Parameter k");
  g1->GetYaxis()->SetTitleOffset(1.5);
  //g1->GetXaxis()->SetLimits(0.0,2*muFitted);
  //g1->GetHistogram()->SetMinimum(0.0);
  //g1->GetHistogram()->SetMaximum(2*kFitted);
  g1->Draw();
  g1->Draw("same L*");
  g2->Draw("same L*");
  g3->Draw("same L*");

  TCanvas * c4 = new TCanvas("c4","c4",400,400,400,400);
  c4->cd();
  c4->SetLeftMargin(0.12);
  c4->SetLogy();
  TGraph *gworst = new TGraph(Nsteps,stepvals,chi2worstvals);
  gworst->SetLineColor(kRed);
  TGraph *gbest = new TGraph(Nsteps,stepvals,chi2bestvals);
  gbest->SetLineColor(kBlue);
  gworst->SetTitle("Best and worst #chi^{2}");
  gworst->GetXaxis()->SetTitleSize(0.04);
  gworst->GetXaxis()->SetLabelSize(0.04);
  gworst->GetYaxis()->SetTitleSize(0.04);
  gworst->GetYaxis()->SetLabelSize(0.04);
  gworst->GetXaxis()->SetTitle("steps");
  //gworst->GetXaxis()->SetTitleOffset(1.2);
  gworst->GetYaxis()->SetTitle("#chi^{2}");
  //gworst->GetYaxis()->SetTitleOffset(1.8);
  gworst->SetMinimum(100.0);
  gworst->SetMaximum(1e4);
  gworst->Draw();
  gworst->Draw("same L");
  gbest->Draw("same L");

  cout << "Generating final histogram" << endl;
  //Show final histogram and create fine histogram to estimate statistics
  c1->cd();
  delete hSumEnergy_100bins;
  hSumEnergy_100bins = (TH1D*)hNtracks_CMSdata->Clone("HSumEnergy100bins");
  hSumEnergy_100bins->SetLineColor(2);
  TH1D* hSumEnergy_nominal = new TH1D ("HSumEnergy","Energy distribution GM",330,0,histomax);
  for (int ipp = 0; ipp<NEvents; ipp++) {
    //For each nucleon collison (Ncol), sample from the NBD to estimate SumEnergy
    SumEnergy = 0;
    GMntuple->GetEntry(ipp);
    numcol = (int)numcolLeaf->GetValue();
    for (int icolcount = 0; icolcount<numcol; icolcount++) {
      SumEnergy += NBD->GetRandom()*xscale;
    }
    hSumEnergy_100bins->Fill(SumEnergy);//+xshift);
    hSumEnergy_nominal->Fill(SumEnergy);//+xshift);
  }
  //hSumEnergy_100bins->Scale(2192968/NEvents);
  //hSumEnergy_100bins->Scale(30938610/NEvents);
  hSumEnergy_100bins->Scale(normval/NEvents);
  cout << "total integral = " << hSumEnergy_100bins->Integral() << endl;
  hSumEnergy_100bins->Draw("same hist");
  TText *kmulabel = new TText();
  kmulabel->SetTextFont(1);
  kmulabel->SetTextColor(1);
  kmulabel->SetTextSize(0.03);
  kmulabel->SetTextAlign(12);
  TString kmutext = Form("k=%.3f, mu=%.3f",k,mu);
  kmulabel->DrawText(150, 5e4, kmutext);


  //Display time elapsed
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
  cout << "The program finished." << endl;

  c1->SaveAs("simplexFit_GlauberHFPbPb_MinBiasDataSumETHF.pdf");
  c2->SaveAs("simplexFitPath.pdf");
  c3->SaveAs("simplexFitProcess.pdf");
  c4->SaveAs("simplexFitChi2.pdf");
}

