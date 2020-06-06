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

void MyFinalFit_MinBiasDataSumETHF() {

  cout << "The program started." << endl;

  clock_t start, end;
  double cpu_time_used;
  start = clock();


//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 1: GENERATE EVENTS WITH GLAUBER MODEL MONTE CARLO


  cout << "Loading Glauber Model Data." << endl;

  TFile *GMdata = new TFile("GlauberModelpPbData10000000.root");
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
  TCanvas * c1 = new TCanvas();
  c1->SetLogy();
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.2);

  //Import CMS data.
  const char *inputFile="pPbActivityHistos_MinBias.root";
  TFile *thefile = new TFile(inputFile);
  TH1D* hNtracks_CMSdata = (TH1D*)thefile->Get("hSumET_HF");
  //hNtracks_CMSdata->Scale(2192968);
  hNtracks_CMSdata->Sumw2();
  //hNtracks_CMSdata->Scale(1.0/2192968);
  hNtracks_CMSdata->Draw("E");
  hNtracks_CMSdata->GetYaxis()->SetTitle("Count");
  hNtracks_CMSdata->GetXaxis()->SetTitle("SumET_HF");
  hNtracks_CMSdata->SetMaximum(1e7);
  hNtracks_CMSdata->SetMinimum(0.3);
  //hNtracks_CMSdata->SetMaximum(2e6);
  //hNtracks_CMSdata->SetMinimum(1e5);
  c1->SetBottomMargin(0.15);
  //c1->SetLeftMargin(0.2);

  cout << "Applying particle production model." << endl;

  //Define negative binomial distribution (NBD) with initial parameters
  double mu = 2.687;
  double k = 0.1248;
  double normval = 34.325*1e6/NEvents;//35 works well with 1000000 events.
  double xscale = 1.2592;
  double xshift = 0;
  TF1 *NBD = new TF1("NBD","TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])) * (([0]/[1])**x)/(([0]/[1]+1)**(x+[1]))",0,28);
  NBD->SetParameters(mu,k);

  //Apply particle production model to get the SumEnergy histogram
  double SumEnergy;
  int numcol;
  double histomax = 330;
  TH1D* hSumEnergy_100bins = new TH1D ("HSumEnergy100bins","Energy distribution in pPb GM",165,0,histomax);
  hSumEnergy_100bins->GetYaxis()->SetTitle("Count");
  hSumEnergy_100bins->SetLineColor(2);

  //Scan over k and mu values to find optimum fit.
  int scansize = 5;
  double kmin = 0.34;
  double kmax = 0.44;
  double keps = (kmax-kmin)/scansize;
  double mumin = 2.6;
  double mumax = 2.8;
  double mueps = (mumax-mumin)/scansize;
  double xscalemin = 3.0;
  double xscalemax = 10;
  double xscaleeps = (xscalemax-xscalemin)/scansize;
  double chi2array[scansize][scansize];
  TH2F* kmuscan = new TH2F ("kmuscan","k-mu scan",scansize,kmin,kmax,scansize,mumin,mumax);
  double ktest[scansize], mutest[scansize];
  //xscale = xscalemin-xscaleeps/2;
 //for (int xinc = 0; xinc<scansize; xinc++) {
  //xscale += xscaleeps;
  k = kmin-keps/2;
  for (int kinc = 0; kinc<scansize; kinc++) {
    k += keps;
    mu = mumin-mueps/2;
    for (int muinc = 0; muinc<scansize; muinc++) {
      mu += mueps;
      NBD->SetParameters(mu,k);
      //delete histogram and recreate it.
      delete hSumEnergy_100bins;
      hSumEnergy_100bins = new TH1D ("HSumEnergy100bins","Energy distribution in pPb GM",165,0,histomax);
      for (int ipp = 0; ipp<NEvents; ipp++) {
        //For each nucleon collison (Ncol), sample from the NBD to estimate SumEnergy
        SumEnergy = 0;
        ntuple->GetEntry(ipp);
        numcol = (int)numcolLeaf->GetValue();
        for (int icolcount = 0; icolcount<numcol; icolcount++) {
          SumEnergy += NBD->GetRandom()*xscale;
        }
        hSumEnergy_100bins->Fill(SumEnergy);
      }

      double chi2 = 0;
      /*int normN = 20;
      double normmin = 0.34;
      double normmax = 0.35;
      double normeps = (normmax-normmin)/normN;
      double norm = normmin-normeps/2;
      double chi2min = (double)NEvents;
      TH1D* hnorm = new TH1D("hnorm","hnorm",10,normmin,normmax);
      for (int inorm = 0; inorm < normN; inorm++) {
        norm += normeps;*/
        hSumEnergy_100bins->Scale(normval);
        chi2 = 0;

        //Calculate chi^2 for current choice of k, <n>.
        double Observed, Expected, diff, error;
        for (int ibin = 10; ibin<165; ibin++) {
          Observed = (double)hNtracks_CMSdata->GetBinContent(ibin);
          Expected = (double)hSumEnergy_100bins->GetBinContent(ibin);
          error = (double)hNtracks_CMSdata->GetBinError(ibin);
          if (error > 0) {
            diff = Observed - Expected;
            chi2 += diff*diff/(error*error);
            hSumEnergy_100bins->SetBinContent(0,ibin);
          }
        }
        cout << "k = " << k << "; mu = " << mu << "; norm = " << normval << "; xscale = " << xscale << "; Chi^2 = " << chi2 << endl;

        /*hnorm->SetBinContent(inorm,chi2);
        hSumEnergy_100bins->Scale(1/norm);
        if (chi2 < chi2min) {
          chi2min = chi2;
          normval = norm;
        }
      }

      TCanvas* cnorm = new TCanvas();
      hnorm->Draw();
      hnorm->Fit("pol2");
      cout << "min chi2 = " << chi2min << " at norm = " << normval << "; " << endl;*/

      chi2array[kinc][muinc] = chi2;
      kmuscan->Fill(k,mu,chi2);

    }//end of mu loop

  }//end of k loop

  //Show results of the scan
  TCanvas * cscan = new TCanvas();
  cscan->cd();
  kmuscan->SetStats(kFALSE);
  kmuscan->Draw("colz");
  kmuscan->GetXaxis()->SetTitle("k");
  kmuscan->GetYaxis()->SetTitle("#mu");
  TText *chilabel = new TText();
  chilabel->SetTextFont(1);
  chilabel->SetTextColor(1);
  chilabel->SetTextSize(0.03);
  chilabel->SetTextAlign(12);
  chilabel->SetTextAngle(0);
  //xscale = xscalemin-xscaleeps/2;
 //for (int xinc = 0; xinc<scansize; xinc++) {
  //xscale += xscaleeps;
  k = kmin-keps/2;
  for (int kinc = 0; kinc<scansize; kinc++) {
    k += keps;
    mu = mumin-mueps/2;
    for (int muinc = 0; muinc<scansize; muinc++) {
    mu += mueps;
      TString chistring = Form("%.1f",chi2array[kinc][muinc]);
      chilabel->DrawText(k, mu, chistring);
    }
  }

  //Find minimum bin
  int muminbin, kminbin, zminbin;
  kmuscan->GetMinimumBin(kminbin, muminbin, zminbin);
  k = kmuscan->GetXaxis()->GetBinCenter(kminbin);
  mu = kmuscan->GetYaxis()->GetBinCenter(muminbin);
  cout << "min is at k = " << k << "; mu = " << mu << "; norm = " << normval << ";" << endl;
  cout << "Recommend next scan in range " << k-keps/2 << " < k < " << k+keps/2 << ", " << mu-mueps/2 << " < mu < " << mu+mueps/2 << ". " << endl;
  NBD->SetParameters(mu,k);

  cout << "Generating final histogram" << endl;
  //Show final histogram and create fine histogram to estimate statistics
  c1->cd();
  delete hSumEnergy_100bins;
  hSumEnergy_100bins = new TH1D ("HSumEnergy100bins","Energy distribution in pPb GM",165,0,histomax);
  hSumEnergy_100bins->SetLineColor(2);
  TH1D* hSumEnergy_nominal = new TH1D ("HSumEnergy","Energy distribution GM",330,0,histomax);
  for (int ipp = 0; ipp<NEvents; ipp++) {
    //For each nucleon collison (Ncol), sample from the NBD to estimate SumEnergy
    SumEnergy = 0;
    ntuple->GetEntry(ipp);
    numcol = (int)numcolLeaf->GetValue();
    for (int icolcount = 0; icolcount<numcol; icolcount++) {
      SumEnergy += NBD->GetRandom()*xscale;
    }
    hSumEnergy_100bins->Fill(SumEnergy+xshift);
    hSumEnergy_nominal->Fill(SumEnergy+xshift);
  }
  //hSumEnergy_100bins->Scale(2192968/NEvents);
  //hSumEnergy_100bins->Scale(30938610/NEvents);
  hSumEnergy_100bins->Scale(normval);
  hSumEnergy_100bins->Draw("same");
  TText *kmulabel = new TText();
  kmulabel->SetTextFont(1);
  kmulabel->SetTextColor(1);
  kmulabel->SetTextSize(0.03);
  kmulabel->SetTextAlign(12);
  TString kmutext = Form("k=%.3f, mu=%.3f",k,mu);
  kmulabel->DrawText(150, 5e4, kmutext);


//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 3: FIND AVERAGE VALUES OF B, NPART, AND NCOL FOR ALL CENTRALITY CLASSES


  cout << "Estimating class statistics." << endl;


  c1->cd();
  //Normalize histogram and shrink error bars
  //TCanvas * c2 = new TCanvas();
  //c2->cd();
  //c2->SetLogy();
  hSumEnergy_nominal->Sumw2();
  hSumEnergy_nominal->GetYaxis()->SetTitle("Fraction of Events/0.05 TeV");
  hSumEnergy_nominal->GetYaxis()->CenterTitle();
  hSumEnergy_nominal->GetXaxis()->SetTitle("#Sigma E_{T} in HF [TeV]");
  hSumEnergy_nominal->GetXaxis()->CenterTitle();
  //c2->SetBottomMargin(0.15);
  //c2->SetLeftMargin(0.2);
  //hSumEnergy_nominal->Draw("E");

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
  const int xarraysize = 100/percentsize - 1;
  double xtrue[xarraysize] = {0};
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
    double y0 = hSumEnergy_nominal->GetBinContent(bmin);
    int b1 = axis->FindBin(xmin+eps);
    double y1 = hSumEnergy_nominal->GetBinContent(b1);
    double Area0 = integral - inttrue;
    double Area1 = inttrue - previous;
    double A = y0-y1;
    double B = -2*(Area1+Area0) - (y0-y1)*(2*xmin+eps);
    double C = 2*xmin*Area1 + 2*(xmin+eps)*Area0 + (y0-y1)*xmin*(xmin+eps);
    xtrue[n-1] = (-B-TMath::Sqrt(B*B-4*A*C))/(2*A);
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
  double xprev = 200;
  for (int i = 0; i<xarraysize; i++) {
    const double xpos = xtrue[i];
    cout << "xtrue[" << i << "] = " << xpos << endl;
    TLine *l = new TLine(xpos,0.3,xpos,1e7);
    l->SetLineColor(2);
    l->SetLineStyle(2);
    l->Draw();
    const double xavg = (xpos+xprev)/2;
    TString percent = Form("%d - %d %%",i*percentsize,(i+1)*percentsize);
    if (i<6) {
      label->DrawText(xavg, 1, percent);
    }
    xprev = xpos;
  }

  TLegend* fitleg = new TLegend(0.76,0.5,0.9,0.7); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(hNtracks_CMSdata,"Data","pel");
  fitleg->AddEntry(hSumEnergy_100bins,"GM","l");
  fitleg->Draw("same");

  cout << "Finished drawing" << endl;

  //Save the Histogram
  c1->SaveAs("GlauberHFpPb_MinBiasDataSumETHF_1.pdf");
  c1->SaveAs("GlauberHFpPb_MinBiasDataSumETHF_1.png");

  cout << "Saved histograms" << endl;

  //Determine average values of Npart, Ncol, and b for each centrality class
  /*double Npartavg[xarraysize+1] = {0};
  double Ncolavg[xarraysize+1] = {0};
  double bavg[xarraysize+1] = {0};
  double b;
  int numpart;
  int eventcounter = 0;
  int classsize[xarraysize+1] = {0};
  for (int eventcounter = 0; eventcounter < NEvents; eventcounter++) {
    int classcounter = -1;
    double test = SumEnergy[eventcounter];
    if (test > xtrue[0]) {
      classcounter = 0;
    }
    for (int j = 1; j<xarraysize; j++) {
      if (xtrue[j-1] > test && test > xtrue[j]) {
        classcounter = j;
      }
    }
    if (classcounter < 0) {
      classcounter = xarraysize;
    }
    ntuple->GetEntry(eventcounter);
    b = (double)bLeaf->GetValue();
    numpart = (int)numpartLeaf->GetValue();
    numcol = (int)numcolLeaf->GetValue();
    Npartavg[classcounter] += numpart;
    Ncolavg[classcounter] += numcol;
    bavg[classcounter] += b;
    classsize[classcounter]++;
  }

  //Diplay average values of Npart, Ncol, and b for each centrality class
  for (int classcounter2 = 0; classcounter2<xarraysize+1; classcounter2++) {
    Npartavg[classcounter2] = (double)Npartavg[classcounter2]/classsize[classcounter2];
    Ncolavg[classcounter2] = (double)Ncolavg[classcounter2]/classsize[classcounter2];
    bavg[classcounter2] = bavg[classcounter2]/classsize[classcounter2];
    cout << "Npartavg[" << classcounter2 << "] = " << Npartavg[classcounter2] << "; ";
    cout << "Ncolavg[" << classcounter2 << "] = " << Ncolavg[classcounter2] << "; ";
    cout << "bavg[" << classcounter2 << "] = " << bavg[classcounter2] << "; ";
    cout << "n = " << classsize[classcounter2] << endl;
  }*/

  //Display time elapsed
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
  cout << "The program finished." << endl;
}

