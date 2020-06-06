//Author: Jared Jay
//This code takes two Ttrees stored in two different data files and outputs a merged Ttree in a new file.

#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"

using namespace std;

void mergeDataTrees() {

  TString f1Name = "../yskimPA1st_OpSign_20177262037_unIdentified.root";
  TString f2Name = "../yskimPA2nd_OpSign_20177262044_unIdentified.root";
  TString outFileName = "../yskimPA1stAnd2nd_OpSign_20177262044_unIdentified.root";

  cout << "OPENING DATA FILES" << endl;
  TFile* f1 = new TFile(f1Name);
  TFile* f2 = new TFile(f2Name);
  TFile* outFile = new TFile(outFileName,"RECREATE");

  TTree* tree1 = (TTree*)f1->Get("myTree");
  TTree* tree2 = (TTree*)f2->Get("myTree");

  cout << "MERGING TREES" << endl;
  TList *list = new TList;
  list->Add(tree1);
  list->Add(tree2);
  TTree *newtree = TTree::MergeTrees(list);
  //newtree->SetName("newtree");

  cout << "WRITING TO FILE" << endl;
  //newtree->Write();
  outFile->Close();

}
