#include "Chain.hh"
#include <iostream>

Chain::Chain(){
}

Chain::Chain(Int_t series, Int_t first, Int_t last, TString dir, 
	     const char* treename){
  fSeries = series;
  fFirst = first;
  fLast = last;
  fDirName = dir;
  fTreeName = (char *)treename;
  CreateChain();
}

Chain::Chain(Int_t series, Int_t *filearray, Int_t len, TString dir, 
	     const char* treename){
  fSeries = series;
  fFileArray = filearray;
  fLen = len;
  fDirName = dir;
  fTreeName = (char *)treename;
  CreateChain();
}

void Chain::CreateChain(){
  TString dname = fDirName + Form("/series%03i", fSeries);
  fChain = new TChain(fTreeName);
  if(fLen > 0){//complete chain with specified array of files
    for(int i=0;i<fLen; ++i){
      fChain->Add(Form("%s/s%03i_f%05i*.root",dname.Data(), fSeries, 
		       fFileArray[i]));
      if(fMaxPulses > 0){
	Long64_t np = fChain->GetEntries();
	if(np > fMaxPulses)break;
      }
    }
  }else{ //complete chain with specified range of files
    int nf = fLast - fFirst + 1;
    for(int i=0;i<nf; ++i){
      TString fn = Form("%s/s%03i_f%05i*.root",dname.Data(), fSeries, fFirst+i);
      int nadd  = fChain->Add(fn.Data());
      fNfiles += nadd;
      if(fVerbose){
	if(nadd)
	  printf("Added s%03i_f%05i*.root\n", fSeries, fFirst+i);
	else
	  printf("s%03i_f%05i*.root not found\n", fSeries, fFirst+i);
      }
      if(fMaxPulses > 0){
	Long64_t np = fChain->GetEntries();
	if(np > fMaxPulses)break;
      }
    }
  }
  printf("%i Files added to chain\n", fNfiles);
  if(fNfiles > 0)
    fIsMade = true;
}

TChain* Chain::GetChain(){
  if(fIsMade){
    return fChain;
  }
  else {
    printf("No TChain created yet. Run CreateChain() first.\n");
    return 0;
  }
}

void Chain::SetDir(TString dir){fDirName = dir;}

void Chain::SetFileArray(Int_t* filearray, int len){
  fLen = len;
  fFileArray = filearray;
}

void Chain::SetFirst(Int_t first){fFirst = first;}

void Chain::SetLast(Int_t last){fLast = last;}

void Chain::SetMaxPulses(Int_t max_pulses){fMaxPulses = max_pulses;}

void Chain::SetSeries(Int_t series){fSeries = series;}

void Chain::SetTreeName(const char* treename){fTreeName = (char *)treename;}

void Chain::SetVerbose(Bool_t verbose){fVerbose = verbose;}


