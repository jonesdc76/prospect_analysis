#include <fstream>
#include <string>
#include <cstdio>
#include <iostream>
#include "stdio.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TPaveText.h"

inline bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}





TChain* ChainIBD( int start = 0 , int end = 1199,
		  TString dir_pref = "Mixed_IBD_with_muBD_and_nBD_",
		  TString dataset = "IBD_with_BG")
{
  //  TCanvas *c = new TCanvas("c","c",0,0,500,600);
  TChain *ch = new TChain("IBDTreePlugin/Tibd");
  int n = 0;
  for(int i=start;i<=end;++i){
    string fname = Form("%s/AD1_%s/%sRun_%i/P2kAnalyzer.root",gSystem->Getenv("DATA_CHALLENGE_OUTDIR"),dataset.Data(),dir_pref.Data(),i);
    if(file_exists(fname.c_str())){
      ch->Add(fname.c_str());
      cout<<fname.c_str()<<endl;
      //cout<<"Run "<<i<<": "<<ch->GetEntries()<<" entries"<<endl;
      ++n;
    }else
      cout<<fname.c_str()<<" not found"<<endl;
  }
  cout<<n<<" runs found"<<endl;

  return ch;
}
