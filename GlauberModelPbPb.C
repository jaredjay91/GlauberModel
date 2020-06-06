#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TH2.h"
#include "time.h"
#include "TNtuple.h"
#include "TFile.h"

void GlauberModelPbPb(){

  cout << "The program started." << endl;

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  //set constant parameters
  const double pi = TMath::Pi();
  const double R = 6.62;
  const double a = 0.546;
  const double max = 2*R;
  const double sigma = 7.0; //70 mb (1 b = 10^-28 m^2 = 100 fm^2)
  const double dtest = sigma/pi;
  const double scale = 5;
  const int NEvents = 1e6;//CMS used 25397.
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
  TNtuple* ntuple = new TNtuple("ntuple","Data from Glauber Model","b:numpart:numcol",NEvents);
  int iEvents = 0;

  //loop until we have enough collisions
  while (iEvents<NEvents) {

    b = bprob->GetRandom();

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
    
    //Move on if there's a collision or else re-do that event
    if (numcol>0) {
      iEvents++;
      ntuple->Fill(b,numpart,numcol);
    }
  }

  //Save data
  TFile outfile ("GlauberModelData1000000.root", "RECREATE");
  ntuple->Write();
  outfile.Close();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
}
