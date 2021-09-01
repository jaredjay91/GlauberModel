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

void getClasses() {

  cout << "The program started." << endl;

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  gStyle->SetOptStat(0);

//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 1: GENERATE EVENTS WITH GLAUBER MODEL MONTE CARLO


  cout << "Loading Glauber Model Data." << endl;

  TFile *GMdata = new TFile("GlauberModelDataPbPb502TeV1000000.root");
  TNtuple* ntuple = (TNtuple*)GMdata->Get("ntuple;1");
  TLeaf *bLeaf = ntuple->GetLeaf("b");
  TLeaf *numpartLeaf = ntuple->GetLeaf("numpart");
  TLeaf *numcolLeaf = ntuple->GetLeaf("numcol");
  const int NEvents = (int)ntuple->GetEntries();
  cout << "NEvents = " << NEvents << endl;
  cout << "Glauber Model loaded." << endl;


//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 2: APPLY PARTICLE PRODUCTION MODEL AND FIT MONTE CARLO TO CMS DATA

  //Create canvas
  TCanvas * c1 = new TCanvas("c1","c1",0,0,500,400);
  c1->SetLogy();
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.2);

  //Import CMS data.
  const char *inputFile = "10-001_Sum_ET_HF_Energy.root";
  TFile *thefile = new TFile(inputFile);
  TH1D* hNtracks_CMSdata = (TH1D*)thefile->Get("hSumEnergy_nominal");
  double Ndata = hNtracks_CMSdata->GetEntries();
  int minBinFit = 5;
  int nHFbins = hNtracks_CMSdata->GetNbinsX();
  cout << "numClasses of bins = " << nHFbins << endl;

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
  //hNtracks_CMSdata->GetYaxis()->SetTitle("Count");
  //hNtracks_CMSdata->GetXaxis()->SetTitle("SumET_HF");
  hNtracks_CMSdata->GetYaxis()->SetTitle("Fraction of Events/0.05 TeV");
  hNtracks_CMSdata->GetYaxis()->CenterTitle();
  hNtracks_CMSdata->GetXaxis()->SetTitle("#Sigma E_{T} in HF [TeV]");
  hNtracks_CMSdata->GetXaxis()->CenterTitle();
  hNtracks_CMSdata->GetYaxis()->SetTitleSize(0.05);
  hNtracks_CMSdata->GetYaxis()->SetLabelSize(0.05);
  hNtracks_CMSdata->GetXaxis()->SetTitleSize(0.05);
  hNtracks_CMSdata->GetXaxis()->SetTitleOffset(1.5);
  hNtracks_CMSdata->GetXaxis()->SetLabelSize(0.05);
  //hNtracks_CMSdata->SetMaximum(1e7);
  //hNtracks_CMSdata->SetMinimum(0.3);
  //hNtracks_CMSdata->SetMaximum(2e6);
  //hNtracks_CMSdata->SetMinimum(1e5);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);

  cout << "total integral = " << hNtracks_CMSdata->Integral() << endl;

  cout << "Applying particle production model." << endl;

  //Define negative binomial distribution (NBD) with initial parameters
  //double mu = 1.44058;
  //double k = 0.782729;
  //double normval = 1.195;
  double mu = 1.32399;
  double k = 0.461721;
  double normval = 1.1975;
  //double xscale = 1.2592;
  double xscale = 1.0/1000;
  //double xshift = 0;
  TF1 *NBD = new TF1("NBD","TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])) * (([0]/[1])**x)/(([0]/[1]+1)**(x+[1]))",0,28);
  NBD->SetParameters(mu,k);

  //Apply particle production model to get the SumEnergy histogram
  double SumEnergy[NEvents];
  int numcol;
  double histomax = 5;

  TH1D* hSumEnergy_100bins = (TH1D*)hNtracks_CMSdata->Clone("HSumEnergy100bins");
  hSumEnergy_100bins->SetLineColor(2);

  cout << "Generating final histogram" << endl;
  //Show final histogram and create fine histogram to estimate statistics
  c1->cd();
  delete hSumEnergy_100bins;
  hSumEnergy_100bins = (TH1D*)hNtracks_CMSdata->Clone("HSumEnergy100bins");
  hSumEnergy_100bins->SetLineColor(2);
  TH1D* hSumEnergy_nominal = new TH1D ("HSumEnergy","Energy distribution GM",1000,0,histomax);
  //hSumEnergy_nominal->SetLineColor(2);
  for (int ipp = 0; ipp<NEvents; ipp++) {
    //For each nucleon collison (Ncol), sample from the NBD to estimate SumEnergy
    SumEnergy[ipp] = 0;
    ntuple->GetEntry(ipp);
    numcol = (int)numcolLeaf->GetValue();
    for (int icolcount = 0; icolcount<numcol; icolcount++) {
      SumEnergy[ipp] += NBD->GetRandom()*xscale;
    }
    hSumEnergy_100bins->Fill(SumEnergy[ipp]);//+xshift);
    hSumEnergy_nominal->Fill(SumEnergy[ipp]);//+xshift);
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


//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 3: FIND AVERAGE VALUES OF B, NPART, AND NCOL FOR ALL CENTRALITY CLASSES

  cout << "Estimating class statistics." << endl;

  c1->cd();
  hSumEnergy_nominal->Sumw2();
  //Integrate in steps of 10% of the area under the curve and keep track of the resulting limits
  TAxis *axis = hSumEnergy_nominal->GetXaxis();
  int bmin = axis->FindBin(0.0);
  int bmax = hSumEnergy_nominal->FindLastBinAbove();
  int nbins = hSumEnergy_nominal->GetSize() - 2;
  int percentsize = 10;
  int n = 1;
  double eps = hSumEnergy_nominal->GetBinWidth(bmax);
  double xmin = histomax;
  double xmax = histomax;
  double totalintegral = hSumEnergy_nominal->Integral(bmin,bmax);
  hSumEnergy_nominal->Scale(1.0/totalintegral);
  hSumEnergy_nominal->SetMaximum(0.03);
  hSumEnergy_nominal->SetMinimum(1e-6);
  const int numClasses = (int)(100/percentsize);
  const int xarraysize = numClasses - 1;
  double xtrue[xarraysize];
  memset( xtrue, 0, xarraysize*sizeof(double) );
  while (xtrue[xarraysize-1]==0) {
    double integral = 0.0;
    double previous = 0.0;
    double inttrue = (double)percentsize*n/100.0;
    while (integral < inttrue) {
      previous = integral;
      xmin -= eps;
      bmin = axis->FindBin(xmin);
      bmax = axis->FindBin(xmax);
      integral = hSumEnergy_nominal->Integral(bmin,bmax);
    }
    //estimate true x-value
    cout << "integral = " << integral << " >= " << inttrue << " at xmin=" << xmin << endl;
    cout << "previous = " << previous << " < " << inttrue << " at xmin=" << xmin+eps << endl;
    double y0 = hSumEnergy_nominal->GetBinContent(bmin);
    int b1 = axis->FindBin(xmin+eps);
    double y1 = hSumEnergy_nominal->GetBinContent(b1);
    double Area0 = integral - inttrue;
    double Area1 = inttrue - previous;
    double A = y0-y1;
    double B = -2*(Area1+Area0) - (y0-y1)*(2*xmin+eps);
    double C = 2*xmin*Area1 + 2*(xmin+eps)*Area0 + (y0-y1)*xmin*(xmin+eps);
    xtrue[n-1] = (-B-TMath::Sqrt(B*B-4*A*C))/(2*A);
    cout << "it should be exactly " << inttrue << " at xtrue=" << xtrue[n-1] << endl;
    n++;
  }
  cout << "Finished while loop" << endl;

  //Draw the limits on the graph
  TText *label = new TText();
  label->SetTextFont(1);
  label->SetTextColor(1);
  label->SetTextSize(0.03);
  label->SetTextAlign(12);
  label->SetTextAngle(90);
  double xprev = 5;
  double ymin = hSumEnergy_100bins->GetMinimum();
  double ymax = hSumEnergy_100bins->GetMaximum();
  double ymid = 1e-4;
  cout << "ymid = " << ymid << endl;
  kmulabel->DrawText(3.5, ymax*0.5, kmutext);
  for (int i = 0; i<xarraysize; i++) {
    const double xpos = xtrue[i];
    cout << "xtrue[" << i << "] = " << xpos << endl;
    TLine *l = new TLine(xpos,ymin,xpos,ymax*0.8);
    l->SetLineColor(2);
    l->SetLineStyle(2);
    l->Draw();
    const double xavg = (xpos+xprev)/2;
    TString percent = Form("%d - %d %%",i*percentsize,(i+1)*percentsize);
    if (i<6) {
      label->DrawText(xavg, ymid, percent);
    }
    xprev = xpos;
  }

  TLegend* fitleg = new TLegend(0.7,0.7,0.89,0.8); fitleg->SetTextSize(16);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(hNtracks_CMSdata,"Data","pel");
  fitleg->AddEntry(hSumEnergy_100bins,"Glauber","l");
  fitleg->Draw("same");

  cout << "Finished drawing" << endl;

  //Save the Histogram
  c1->SaveAs(Form("GlauberClasses_PbPbMinBiasDataSumETHF_%iperc.pdf",percentsize));
  c1->SaveAs(Form("GlauberClasses_PbPbMinBiasDataSumETHF_%iperc.png",percentsize));
  cout << "Saved histograms" << endl;

  //Determine average values of Npart, Ncol, and b for each centrality class
  TH1D* hNpartTot = new TH1D("hNpartTot","Number of participants",90,0,450);
  TH1D* hNcolTot = new TH1D("hNcolTot","Number of Collisions",100,0,2500);
  TH1D* hbTot = new TH1D("hbTot","Impact Parameter",100,0,20);
  TH1D* hNpart[numClasses];
  TH1D* hNcol[numClasses];
  TH1D* hb[numClasses];
  for (int iClass = 0; iClass<numClasses; iClass++) {
    hNpart[iClass] = new TH1D(Form("hNpart[%i]",iClass),"hNpart",90,0,450);
    hNcol[iClass] = new TH1D(Form("hNcol[%i]",iClass),"hNcol",100,0,2500);
    hb[iClass] = new TH1D(Form("hb[%i]",iClass),"hb",100,0,20);
  }
  double Npartavg[numClasses];
  double Ncolavg[numClasses];
  double bavg[numClasses];
  memset( Npartavg, 0, numClasses*sizeof(double) );
  memset( Ncolavg, 0, numClasses*sizeof(double) );
  memset( bavg, 0, numClasses*sizeof(double) );
  double b;
  int numpart;
  int iEvt = 0;
  int classSize[numClasses];
  memset( classSize, 0, numClasses*sizeof(double) );
  for (int iEvt = 0; iEvt < NEvents; iEvt++) {
    int iClass = -1;
    double test = SumEnergy[iEvt];
    if (test > xtrue[0]) {
      iClass = 0;
    }
    for (int j = 1; j<xarraysize; j++) {
      if (xtrue[j-1] > test && test > xtrue[j]) {
        iClass = j;
      }
    }
    if (iClass < 0) {
      iClass = xarraysize;
    }
    ntuple->GetEntry(iEvt);
    b = (double)bLeaf->GetValue();
    numpart = (int)numpartLeaf->GetValue();
    numcol = (int)numcolLeaf->GetValue();
    hNpart[iClass]->Fill(numpart);
    hNcol[iClass]->Fill(numcol);
    hb[iClass]->Fill(b);
    hNpartTot->Fill(numpart);
    hNcolTot->Fill(numcol);
    hbTot->Fill(b);
    Npartavg[iClass] += numpart;
    Ncolavg[iClass] += numcol;
    bavg[iClass] += b;
    classSize[iClass]++;
  }

  //Diplay average values of Npart, Ncol, and b for each centrality class
  for (int iClass = 0; iClass<numClasses; iClass++) {
    Npartavg[iClass] = (double)Npartavg[iClass]/classSize[iClass];
    Ncolavg[iClass] = (double)Ncolavg[iClass]/classSize[iClass];
    bavg[iClass] = bavg[iClass]/classSize[iClass];
    cout << "Npartavg[" << iClass << "] = " << Npartavg[iClass] << "; ";
    cout << "Ncolavg[" << iClass << "] = " << Ncolavg[iClass] << "; ";
    cout << "bavg[" << iClass << "] = " << bavg[iClass] << "; ";
    cout << "n = " << classSize[iClass] << endl;
  }

  //LaTeX format
  cout << endl << "Centrality & $\\langle N_{part}\\rangle$ & $\\langle N_{coll}\\rangle$ & b \\\\" << endl;
  for (int iClass = 0; iClass<numClasses; iClass++) {
    cout << Form("%d - %d",iClass*percentsize,(iClass+1)*percentsize) << "\\% & " << Npartavg[iClass] << " & " << Ncolavg[iClass] << " & " << bavg[iClass] << " \\\\" << endl;
  }
  //Compare to https://twiki.cern.ch/twiki/bin/viewauth/CMS/Glauber5TeVPbPbNewParameters#Npart

  //make a fancy looking plot.
  float textSize = 0.05;
  TCanvas* c2 =  new TCanvas("c2","c2",0,0,1200,400);
  c2->Divide(3,1);
  c2->cd(1);
  gPad->SetLogy();
  gPad->SetBottomMargin(0.11);
  hNpartTot->GetXaxis()->SetTitle("N_{part}");
  hNpartTot->GetXaxis()->SetLabelSize(textSize);
  hNpartTot->GetXaxis()->SetTitleSize(textSize);
  hNpartTot->GetYaxis()->SetLabelSize(textSize);
  hNpartTot->GetYaxis()->SetTitleSize(textSize);
  hNpartTot->SetMaximum(0.2*NEvents);
  hNpartTot->SetMinimum(0.6);
  hNpartTot->Draw();
  c2->cd(2);
  gPad->SetLogy();
  gPad->SetBottomMargin(0.11);
  hNcolTot->GetXaxis()->SetTitle("N_{coll}");
  hNcolTot->GetXaxis()->SetLabelSize(textSize);
  hNcolTot->GetXaxis()->SetTitleSize(textSize);
  hNcolTot->GetYaxis()->SetLabelSize(textSize);
  hNcolTot->GetYaxis()->SetTitleSize(textSize);
  hNcolTot->SetMaximum(0.4*NEvents);
  hNcolTot->SetMinimum(0.6);
  hNcolTot->Draw();
  c2->cd(3);
  TGaxis::SetMaxDigits(3);
  gPad->SetBottomMargin(0.11);
  hbTot->GetXaxis()->SetTitle("b");
  hbTot->GetXaxis()->SetLabelSize(textSize);
  hbTot->GetXaxis()->SetTitleSize(textSize);
  hbTot->GetYaxis()->SetLabelSize(textSize);
  hbTot->GetYaxis()->SetTitleSize(textSize);
  hbTot->Draw();
  TLegend* leg = new TLegend(0.5, 0.65, 0.67, 0.89);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  TLegend* leg2 = new TLegend(0.7, 0.65, 0.89, 0.89);
  leg2->SetTextSize(0.04);
  leg2->SetBorderSize(0);
  for (int iClass = 0; iClass<numClasses; iClass++) {
    int iColor = iClass+2;
    if (iColor > 9) iColor = iColor+18;
    c2->cd(1);
    hNpart[iClass]->SetLineColor(iColor);
    hNpart[iClass]->Draw("same");
    c2->cd(2);
    hNcol[iClass]->SetLineColor(iColor);
    hNcol[iClass]->Draw("same");
    if (iClass>4) leg2->AddEntry(hNcol[iClass],Form("%i-%i%%",iClass*10,(iClass+1)*10),"l");
    else leg->AddEntry(hNcol[iClass],Form("%i-%i%%",iClass*10,(iClass+1)*10),"l");
    c2->cd(3);
    hb[iClass]->SetLineColor(iColor);
    hb[iClass]->Draw("same");
  }


  c2->cd(2);
  leg->Draw("same");
  leg2->Draw("same");
  c2->SaveAs(Form("NpartNcolb_centralityClasses_%iperc.png",percentsize));
  c2->SaveAs(Form("NpartNcolb_centralityClasses_%iperc.pdf",percentsize));

  //Display time elapsed
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
  cout << "The program finished." << endl;

}

