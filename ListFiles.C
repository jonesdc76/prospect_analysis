#include <iostream>
#include <list>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;
struct Pulse_t{
  Long64_t evt;
  Double_t t;
  Int_t seg;
  Float_t E;
  Float_t dt;
  Float_t y;
  Float_t PE[2];
  Float_t PSD;
  Int_t PID;
};


void ListFiles(TString fname, TString series = "007", TString data_set = "P50D",
	     TString data_release = "Phys_20161206", 
	     TString conv_data_dir = "/home/jonesdc/prospect/data/Analyzed/") 
{
  FILE *f = fopen("/home/jonesdc/prospect/macros/runlist.dat","a+");
  if(f==0){
    cout<<"File not opened. Exiting.\n";
    return;
  }
 //Open the Root tree
  TString fn = conv_data_dir + "/" + data_release + "/" + data_set + "/" + 
    "series" + series + "/" + fname;
  //Open file with BiPo tree
  fn = gSystem->Getenv("P2X_ANALYZED");
  fn += "/" + data_release + "/" + data_set + "/" + 
    "series" + series; 
  fn += "/biPoEvents_" + fname;
  TFile *file = TFile::Open(fn.Data());
  if(file == 0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return;
  }
  TTree *tr = (TTree*)file->Get("BiPo");
  TCut cut = TCut("");
  double nfar, nprompt;
  if(tr->GetEntries()>0){
    tr->Draw("nFar>>h(100)",cut,"goff");
    TH1D *hf = (TH1D*)gDirectory->Get("h");
    nfar = hf->GetMean();
    tr->Draw("nPrompt>>h1(100)",cut,"goff");
    TH1D *h1 = (TH1D*)gDirectory->Get("h1");
    nprompt = h1->GetMean();
    cout<<fname.Data()<<" "<<nfar<<" "<<nprompt<<" "<<nfar/nprompt<<endl;
    fprintf(f,"%s  %0.5f  %0.5f  %0.5f\n", fname.Data(), nfar, nprompt, nfar/nprompt);
  }
  fclose(f);
  return;

}
