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

void MyFinalFit() {

  cout << "The program started." << endl;

  clock_t start, end;
  double cpu_time_used;
  start = clock();


//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 1: GENERATE EVENTS WITH GLAUBER MODEL MONTE CARLO


  cout << "Loading Glauber Model Data." << endl;

  TFile *GMdata = new TFile("GlauberModelData1000000.root");
  TNtuple* ntuple = (TNtuple*)GMdata->Get("ntuple;1");
  TLeaf *bLeaf = ntuple->GetLeaf("b");
  TLeaf *numpartLeaf = ntuple->GetLeaf("numpart");
  TLeaf *numcolLeaf = ntuple->GetLeaf("numcol");
  const int NEvents = (int)ntuple->GetEntries();

  cout << "NEvents = " << NEvents << endl;

  cout << "Glauber Model loaded." << endl;


//---------------------------------------------------------------------------------------------------------------------------------------------
//PART 2: APPLY PARTICLE PRODUCTION MODEL AND FIT MONTE CARLO TO CMS DATA


  cout << "Applying particle production model." << endl;

  //Define negative binomial distribution (NBD) with initial parameters
  double mu = 2.53;
  double k = 2.00;
  TF1 *NBD = new TF1("NBD","TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1])) * (([0]/[1])**x)/(([0]/[1]+1)**(x+[1]))",0,28);
  NBD->SetParameters(mu,k);

  //Create canvas
  TCanvas * c1 = new TCanvas();
  c1->SetLogy();

  //Import CMS data.
  const char *inputFile="10-001_Sum_ET_HF_Energy.root";
  TFile *thefile = new TFile(inputFile);
  TH1D* hSumEnergy_CMSdata = (TH1D*)thefile->Get("hSumEnergy_nominal;1");
  hSumEnergy_CMSdata->Scale(25397);
  hSumEnergy_CMSdata->Sumw2();
  hSumEnergy_CMSdata->Scale(1.0/25397);
  hSumEnergy_CMSdata->Draw("E");
  hSumEnergy_CMSdata->GetYaxis()->SetTitle("Fraction of Events/0.05 TeV");
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.2);

  //Apply particle production model to get the SumEnergy histogram
  double SumEnergy[NEvents];
  int scansize = 5;
  int numcol;
  double kmin = 0.88;
  double kmax = 0.92;
  double keps = (kmax-kmin)/scansize;
  double mumin = 1.68;
  double mumax = 1.72;
  double mueps = (mumax-mumin)/scansize;
  double chi2array[scansize][scansize];
  TH2F* kmuscan = new TH2F ("kmuscan","k-mu scan",scansize,kmin,kmax,scansize,mumin,mumax);
  double ktest[scansize], mutest[scansize];
  const double scale = 5;//CMS mean = 0.9701
  const double scaleden = 5544;//CMS mean = 0.9701. My data mean = 0.7684*7000/5=1075.76. 7200 works well.
  TH1D* hSumEnergy_100bins = new TH1D ("HSumEnergy100bins","Energy distribution with 100 bins",100,0,scale);
  double normval = 1.0;

  //Scan over k and mu values to find optimum fit.
  k = kmin-keps/2;
  for (int kinc = 0; kinc<scansize; kinc++) {
    k += keps;
    mu = mumin-mueps/2;
    for (int muinc = 0; muinc<scansize; muinc++) {
      mu += mueps;
      NBD->SetParameters(mu,k);
      for (int ipp = 0; ipp<NEvents; ipp++) {
        //For each nucleon collison (Ncol), sample from the NBD to estimate SumEnergy
        SumEnergy[ipp] = 0;
        ntuple->GetEntry(ipp);
        numcol = (int)numcolLeaf->GetValue();
        for (int icolcount = 0; icolcount<numcol; icolcount++) {
          SumEnergy[ipp] += NBD->GetRandom()*scale/scaleden;
        }
        hSumEnergy_100bins->Fill(SumEnergy[ipp]);
      }

      hSumEnergy_100bins->Scale(1.0/NEvents);
      int normN = 1;
      double normmin = 1.225;
      double normmax = 1.230;
      double normeps = (normmax-normmin)/normN;
      double norm = normmin-normeps/2;
      double chi2 = 0;
      double chi2min = 1000.0;
      for (int inorm = 0; inorm < normN; inorm++) {
        norm += normeps;
        hSumEnergy_100bins->Scale(norm);

        //Calculate chi^2 for current choice of k, <n>.
        chi2 = 0;
        double Observed, Expected, diff, error;
        for (int ibin = 10; ibin<100; ibin++) {
          Observed = (double)hSumEnergy_CMSdata->GetBinContent(ibin);
          Expected = (double)hSumEnergy_100bins->GetBinContent(ibin);
          error = (double)hSumEnergy_CMSdata->GetBinError(ibin);
          if (error > 0) {
            diff = Observed - Expected;
            chi2 += diff*diff/(error*error);
            hSumEnergy_100bins->SetBinContent(0,ibin);
          }
        }
        cout << "k = " << k << "; mu = " << mu << "; norm = " << norm << "; Chi^2 = " << chi2 << endl;
        hSumEnergy_100bins->Scale(1/norm);
        if (chi2 < chi2min) {
          chi2min = chi2;
          normval = norm;
        }
      }

      cout << "min chi2 = " << chi2min << " at norm = " << normval << "; " << endl;

      chi2array[kinc][muinc] = chi2;
      kmuscan->Fill(k,mu,chi2);

    }//end of mu loop

  }//end of k loop

  //Find minimum bin
  int muminbin, kminbin, zminbin;
  kmuscan->GetMinimumBin(muminbin, kminbin, zminbin);
  k = kmuscan->GetXaxis()->GetBinCenter(muminbin);
  mu = kmuscan->GetYaxis()->GetBinCenter(kminbin);
  cout << "min is at k = " << k << " and mu = " << mu << endl;
  cout << "Recommend next scan in range " << k-keps/2 << " < k < " << k+keps/2 << ", " << mu-mueps/2 << " < mu < " << mu+mueps/2 << ". " << endl;
  NBD->SetParameters(mu,k);

  //Show final histogram and create fine histogram to estimate statistics
  TH1D* hSumEnergy_nominal = new TH1D ("HSumEnergy","Energy distribution",700,0,scale);
  for (int ipp = 0; ipp<NEvents; ipp++) {
    //For each nucleon collison (Ncol), sample from the NBD to estimate SumEnergy
    SumEnergy[ipp] = 0;
    ntuple->GetEntry(ipp);
    numcol = (int)numcolLeaf->GetValue();
    for (int icolcount = 0; icolcount<numcol; icolcount++) {
      SumEnergy[ipp] += NBD->GetRandom()*scale/scaleden;
    }
    hSumEnergy_100bins->Fill(SumEnergy[ipp]);
    hSumEnergy_nominal->Fill(SumEnergy[ipp]);
  }
  hSumEnergy_100bins->Scale(1.0/NEvents);
  hSumEnergy_100bins->Draw("same");

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
  k = kmin-keps/2;
  for (int kinc = 0; kinc<scansize; kinc++) {
    k += keps;
    mu = mumin-mueps/2;
    for (int muinc = 0; muinc<scansize; muinc++) {
    mu += mueps;
      TString chistring = Form("%f",chi2array[kinc][muinc]);
      chilabel->DrawText(k, mu, chistring);
    }
  }


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
  //int bmax = axis->FindBin(scale);
  int nbins = hSumEnergy_nominal->GetSize() - 2;
  int percentsize = 10;
  int n = 1;
  double eps = hSumEnergy_nominal->GetBinWidth(bmax);
  double xmin = scale;
  double xmax = scale;
  double totalintegral = hSumEnergy_nominal->Integral(bmin,bmax);
  hSumEnergy_nominal->Scale(1.0/totalintegral);
  hSumEnergy_nominal->SetMaximum(0.3);
  hSumEnergy_nominal->SetMinimum(1e-5);
  const int xarraysize = 100/percentsize - 1;
  double xtrue[xarraysize];
  while (xmin > eps) {
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

  hSumEnergy_100bins->Scale(normval);
  //Draw the limits on the graph
  TText *label = new TText();
  label->SetTextFont(1);
  label->SetTextColor(1);
  label->SetTextSize(0.03);
  label->SetTextAlign(12);
  label->SetTextAngle(90);
  double xprev = scale;
  for (int i = 0; i<xarraysize; i++) {
    const double xpos = xtrue[i];
    cout << "xtrue[" << i << "] = " << xpos << endl;
    TLine *l = new TLine(xpos,1e-5,xpos,0.3);
    l->SetLineColor(2);
    l->SetLineStyle(2);
    l->Draw();
    const double xavg = (xpos+xprev)/2;
    TString percent = Form("%d - %d %%",i*percentsize,(i+1)*percentsize);
    if (i<6) {
      label->DrawText(xavg, 1e-4, percent);
    }
    xprev = xpos;
  }

  //Save the Histogram
  c1->SaveAs("HNtupleCentrality.pdf");

  //Determine average values of Npart, Ncol, and b for each centrality class
  double Npartavg[xarraysize+1] = {0};
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
  }

  //Display time elapsed
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
  cout << "The program finished." << endl;
}

