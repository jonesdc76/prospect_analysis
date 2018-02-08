//Class to create a TChain from TFiles from a particular series.
// This chain can be created by providing the first and last file numbers
//to be added and all numbers in between will be assumed, OR by providing 
//an array of file numbers to be added. Alternately, if fLast is 
//not set, fNfiles will be added to the chain starting with fFirst.
//File numbers are from 0 to 99999 and always follow directly after "_f" 
//in the file name. 0 padding is automatic. If desired, the TChain can be
// limited by setting max_pulses.
//Author: Don Jones
//Date: 01/2018
#include "TChain.h"
#include "TString.h"

class Chain{

  private:
  Bool_t fIsMade = false;//set to true when Chain successfully created
  Bool_t fVerbose = 0;//write files added to terminal
  Int_t fSeries = 0;//Series number
  Int_t fFirst = 0;//File number of first file to use in series  
  Int_t fLast = 0; //File number of last file to use in series 
  Long64_t fMaxPulses = 0;//Stop adding files after max_pulses added
  Int_t *fFileArray;//Array of file numbers of interest
  Int_t fLen = 0;//length of filearray
  Int_t fNfiles = 0; //number of files added to TChain
  char* fTreeName;//Name of TTree to load
  TChain *fChain;
  TString fDirName;//Directory where series folders are found

  public:
  Chain();
  Chain(Int_t series, Int_t first, TString dir, char* treename);
  Chain(Int_t series, Int_t first, Int_t last, TString dir, const char* treename);
  Chain(Int_t series, Int_t *filearray, Int_t len, TString dir, const char* treename);
  void CreateChain();
  TChain* GetChain();
  Int_t GetNfiles(){return fNfiles;};
  void SetDir(TString dir);//Directory where files are found
  void SetFileArray(Int_t *filearray, int len);//Array of file numbers to add to chain
  void SetFirst(Int_t first);//File number of first file to use in series  
  void SetLast(Int_t last);//File number of last file to use in series 
  void SetMaxPulses(Int_t n);//Stop adding files after n pulses added
  void SetSeries(Int_t series);//Series number
  void SetTreeName(const char* treename);//Name of TTree
  void SetVerbose(Bool_t verbose);
};

