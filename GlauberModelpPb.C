#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "time.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"

void GlauberModelpPb(){

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
  const int NEvents = 1e8;
  const int numNucleons = 208;

  //create probability distributions
  TF1 *prob = new TF1("prob","(x**2)/(1+TMath::Exp((x-[0])/[1]))",0,max);
  prob->SetParameters(R,a);
  TF1 *thprob = new TF1("thprob","sin(x)",0,pi);
  TF1 *bprob = new TF1("bprob","x",0,1.9*R);

  //initialize variables
  double x1[numNucleons], y1[numNucleons];
  double x2;
  double r, theta, phi, dx, dy, b;
  TRandom3 rnd3(0);
  int numpart, numcol;
  TNtuple* ntuple = new TNtuple("ntuple","Data from Glauber Model","b:numpart:numcol",NEvents);
  int iEvents = 0;

  //loop until we have enough collisions
  while (iEvents<NEvents) {

    b = bprob->GetRandom();
    x2 = b;

    //populate nucleus with nucleons
    for (int i = 0; i<numNucleons; i++) {
      r = prob->GetRandom();
      theta = thprob->GetRandom();
      phi = 2*pi*rnd3.Rndm();
      x1[i] = r*TMath::Sin(theta)*TMath::Cos(phi);
      y1[i] = r*TMath::Sin(theta)*TMath::Sin(phi);
    }

    //search for collisions and count participants
    numcol = 0;
    numpart = 0;
    int participants1[numNucleons] = {0};
    for (int i = 0; i<numNucleons; i++) {
      dx = x1[i]-x2;
      dy = y1[i];
      if (dx*dx + dy*dy < dtest) {
        numcol++;
        numpart++;
        participants1[i] = 1;
      }
    }
    
    //Move on if there's a collision or else re-do that event
    if (numcol>0) {
      numpart++;//because the proton is a participant too.
      iEvents++;
      ntuple->Fill(b,numpart,numcol);
    }
  }

  //Save data
  TString outfileName = Form("GlauberModelpPbData%i.root",NEvents);
  TFile outfile (outfileName, "RECREATE");
  ntuple->Write();
  outfile.Close();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "time elapsed: " << cpu_time_used << " seconds." << endl;
}
