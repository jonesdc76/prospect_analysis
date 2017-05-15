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

//number of cells
const int kNCells = 2;

//half life of Po214 in ns
const double kPo214HalfLife = 164000;

//earliest time (ns) to look for correlated partner 
const double kStart = 10000; 

//latest time (ns) to look for correlated partner
const double kEnd = kStart + 3 * kPo214HalfLife; 

//start time of far window for finding accidental correlated partners
const double kFarStart = kEnd + 10 * kPo214HalfLife; 

//factor by which far window is longer than near
const double kScale = 4.0; 

//latest time to look for correlated partner
const double kFarEnd = kFarStart + kScale * (kEnd - kStart); 



/// stackable! hit PID classifications that can be used in s_PhysPulse
enum HitTypeID {
  //electron/gamma/muon ionization region (might stack with neutron capture)
  IONI_HIT   = 1 << 0,   
     
  //neutron capture (on 6Li) region (might stack with gamma or recoil)
  NCAPT_HIT  = 1 << 1, 
  
  //nucleus recoil (fast neutron) clearly separated from gamma (might stack with neutron capture)
  RECOIL_HIT = 1 << 2, 
  
  //veto detector hit
  VETO_HIT   = 1 << 3,    
  
  //dead volume or non-detected hit (simulation only)    
  DEAD_HIT   = 1 << 4,        

  //unclassifiable event
  CRAZY_HIT  = 1 << 5,        

  //BiPo alpha region hit
  BIPO_HIT   = 1 << 6       
};


struct Pulse_t{
  Long64_t evt;
  Int_t seg;
  Float_t E;
  Double_t t;
  Float_t dt;
  Float_t PE[2];
  Float_t PSD;
  Int_t PID;
};

struct BiPo_t{
  int size; //if prompt filled size=1, prompt and delayed filled>=2,
  int mult_d;//multiplicity of delayed pulses
  int mult_f;//multiplicity of far pulses
  Pulse_t prompt;
  vector<Pulse_t> delayed;
  vector<Pulse_t> far;
};



void test(TString fname, TString series = "007", TString data_set = "P50D",
	     TString data_release = "Phys_20161206", 
	     TString conv_data_dir = "/home/jonesdc/prospect/data/Analyzed/") 
{
  FILE *f = fopen("/home/jonesdc/prospect/macros/FartoNear.dat","a+");
  if(f==0){
    cout<<"File not opened. Exiting.\n";
    return;
  }
  Pulse_t nullevt;
  nullevt.evt = -99999;
  nullevt.seg = -99999;
  nullevt.E = -99999;
  nullevt.t = -99999;
  nullevt.dt = -99999;
  nullevt.PE[0] = -99999;
  nullevt.PE[1] = -99999;
  nullevt.PSD = -99999;
  nullevt.PID = -99999;
  //Open the Root tree
  TString fn = conv_data_dir + "/" + data_release + "/" + data_set + "/" + 
    "series" + series + "/" + fname;
  //Open file with BiPo tree
  fn = gSystem->Getenv("P2X_ANALYZED");
  fn += "/" + data_release + "/" + data_set + "/" + 
    "series" + series; 
  fn += "/BiPoEvents_" + fname;
  TFile *file = TFile::Open(fn.Data());
  if(file == 0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return;
  }
  TTree *tr = (TTree*)file->Get("BiPo");
  TCut cut = TCut("");
  tr->Draw("nFar>>h(100)",cut,"goff");
  TH1D *hf = (TH1D*)gDirectory->Get("h");
  double nfar = hf->GetMean();
  tr->Draw("nDel>>h1(100)",cut,"goff");
  TH1D *hd = (TH1D*)gDirectory->Get("h1");
  double ndel = hd->GetMean();

  cout<<fname.Data()<<" "<<nfar<<" "<<ndel<<" "<<nfar/ndel<<endl;
  fprintf(f,"%s:  %0.5f  %0.5f  %0.5f\n", fname.Data(), nfar, ndel, nfar/ndel);
  fclose(f);
  return;

}
