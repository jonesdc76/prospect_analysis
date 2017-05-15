#include <iostream>
#include <list>
#include <vector>

#include "TFile.h"
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



int BiPoTree_new(TString fname, TString series = "007", TString data_set = "P50D",
	     TString data_release = "Phys_20161206", 
	     TString conv_data_dir = "/home/prospect-collab/data/converted_data") 
{
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
  TFile *file = TFile::Open(fn.Data());
  if(file==0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }
  TTree *tree = (TTree*)file->Get("PhysPulse");

  if(tree == 0){
    cout<<"Cannot find requested ROOT tree in "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }

  Pulse_t pulse; 
  tree->SetBranchAddress("evt", &pulse.evt);
  tree->SetBranchAddress("seg", &pulse.seg);
  tree->SetBranchAddress("E", &pulse.E);  
  tree->SetBranchAddress("t", &pulse.t);
  tree->SetBranchAddress("dt", &pulse.dt);
  tree->SetBranchAddress("PE", &pulse.PE);
  tree->SetBranchAddress("PSD", &pulse.PSD);
  tree->SetBranchAddress("PID", &pulse.PID);

  //Create new ROOT tree to hold BiPo Events
  fn = gSystem->Getenv("P2X_ANALYZED");
  fn += "/" + data_release + "/" + data_set + "/" + 
    "series" + series; 
  //if the directory doesn't exist create it
  TString mkdircmd = "mkdir -p " + fn;
  system(mkdircmd);
  fn += "/BiPoEvents_" + fname;
  TFile *newfile = TFile::Open(fn.Data(), "RECREATE");
  if(newfile == 0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }

  Int_t nFar = 0, nDel = 0;
  Pulse_t prompt;
  Pulse_t Delayed[20], Far[20];
  TTree *BiPoTree = new TTree("BiPo", "Tree containing Bi-Po candidate events");
  BiPoTree->Branch("nFar", &nFar, "nFar/I");
  BiPoTree->Branch("nDel", &nDel, "nDel/I");
  BiPoTree->Branch("prompt", &prompt, "evt/L:seg/I:E/F:t/D:dt/F:PE[2]/F:y/F:PSD/F:PID/I");
  BiPoTree->Branch("delayed", &Delayed, "delayed[nDel]");
  BiPoTree->Branch("far", &Far, "far[nFar]");


  BiPo_t bipo[kNCells];
  for(int i=0; i<kNCells;++i){
    bipo[i].size = 0;
    bipo[i].mult_d = 0;
    bipo[i].mult_f = 0;
  }
  vector<BiPo_t> vBiPo;
  int N = 0;
  tree->GetEntry(0);
  long s_t = pulse.t;
  for(int i=0;i<tree->GetEntries();++i){
    tree->GetEntry(i);
    int seg = pulse.seg;
    if(i%1000==0 && 0)
      cout<<"<<"<<pulse.evt<<" "<<pulse.seg<<" "<<pulse.t<<" "<<pulse.dt
	  <<" "<<pulse.E<<" "<<pulse.PE[0]<<" "<<pulse.PE[1]<<" "
	  <<pulse.PID<<" "<<pulse.PSD<<endl;

    //look for beta-like particles as prompt candidates
    if(pulse.PID & IONI_HIT && bipo[seg].size == 0){
      bipo[seg].prompt = pulse;

      if(0)cout<<"Entry "<<i<<" is beta candidate  in segment "<<seg<<endl;
      ++bipo[seg].size;
      continue;
    }

    //if prompt candidate exists look for delayed candidate(s)
    if(bipo[seg].size>0 && pulse.PID & RECOIL_HIT){//nuclear recoil only
      double d_t = pulse.t - bipo[seg].prompt.t;
      if(d_t > kStart && d_t < kEnd){
	bipo[seg].delayed.push_back(pulse);
	++bipo[seg].size;
	++bipo[seg].mult_d;
	cout<<series<<"/"<<fname<<"  BiPo candidate: "
	    <<bipo[seg].prompt.evt<<" "<<bipo[seg].delayed.back().evt<<" "
	    <<bipo[seg].prompt.seg<<" "<<bipo[seg].delayed.back().seg<<" "
	    <<1e-6*(bipo[seg].prompt.t-s_t)<<" "
	    <<1e-6*bipo[seg].delayed.back().t-s_t<<" "
	    <<bipo[seg].prompt.E<<" "<<bipo[seg].delayed.back().E<<" "
	    <<bipo[seg].prompt.PID<<" "<<bipo[seg].delayed.back().PID<<" "
	    <<bipo[seg].prompt.PSD<<" "<<bipo[seg].delayed.back().PSD<<endl;
      }
      if(d_t > kEnd){//If past correlation window
	if(bipo[seg].size > 1){//if delayed exists
	  if(0)cout<<"Delayed exists: "<<N<<" "<<bipo[seg].mult_d<<endl;
	  //1. check for candidates in time-displaced window
	  int n = 0;
	  while(pulse.t - bipo[seg].prompt.t < kFarEnd){
	    ++N;
	    tree->GetEntry(i+n);
	    if(pulse.t - bipo[seg].prompt.t > kFarStart && pulse.seg == seg){
	      bipo[seg].far.push_back(pulse);
	      ++bipo[seg].mult_f;
	    }
	    ++n;
	  }
	  //2. Transfer from temporary containers to Branch address pointers
	  prompt = bipo[seg].prompt;
	  for(n=0; n < (int)bipo[seg].mult_d;++n)
	    vDelayed.push_back(bipo[seg].delayed[n]);
	  if(bipo[seg].mult_f > 0){
	    for(n=0; n < (int)bipo[seg].mult_f;++n)
	      vFar.push_back(bipo[seg].far[n]);
	  }else
	    vFar.push_back(nullevt);
	  //nDel = bipo[seg].mult_d;
	  //nFar = bipo[seg].mult_f;
	  //cout<<"Del "<<nDel<<"     Far "<<nFar<<endl;
	  //3. Fill the vector
	  //BiPoTree->Fill();
	  vBiPo.push_back(bipo[seg]);
	  if((int)vBiPo.size()>100){
	    for(int i=0;i<(int)vBiPo.size();++i){
	      vFar.clear();
	      vDelayed.clear();
	      nDel = vBiPo[i].mult_d;
	      nFar = vBiPo[i].mult_f;
	      prompt = vBiPo[i].prompt;
	      for(int m=0;m<vBiPo[i].mult_d;++m)
		vDelayed.push_back(vBiPo[i].delayed[m]);
	      for(int m=0;m<vBiPo[i].mult_f;++m)
		vFar.push_back(vBiPo[i].far[m]);
	      //BiPoTree->Fill();
	    }
	    vBiPo.clear();
	  }
	}
	//Reset temporary containers
	bipo[seg].delayed.resize(0);
	bipo[seg].far.resize(0);
	bipo[seg].size = 0;
	bipo[seg].mult_d = 0;
	bipo[seg].mult_f = 0;
	nDel = 0;
	nFar = 0;
      }
    }
  }

  for(int i=0;i<(int)vBiPo.size();++i){
    vFar.clear();
    vDelayed.clear();
    nDel = vBiPo[i].mult_d;
    nFar = vBiPo[i].mult_f;
    prompt = vBiPo[i].prompt;
    for(int m=0;m<vBiPo[i].mult_d;++m)
      vDelayed.push_back(vBiPo[i].delayed[m]);
    for(int m=0;m<vBiPo[i].mult_f;++m)
      vFar.push_back(vBiPo[i].far[m]);
    BiPoTree->Fill();
  }

  BiPoTree->Write("", TObject::kOverwrite);
  newfile->Close();
  file->Close();
  return 0;

}
