#include <ctime>
#include "Chain.cc"
#include "TF1.h"
#include "TFitResult.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <vector>
#include <string.h>

const int kNmaxSamp = 600;
//Args:
//*****************************************************************
//   max_entries   stop after reading this many entries (0 = read all entries)
//   single_pulse  view single pulse trace at a time
//   series        series number
//   first         first file number in series to add to TChain
//   last          last is last file number in series to add to TChain
//   treename      name of tree of which to make TChain
//   dir           directory where series folders are located

int PlotWaveform(Long64_t max_entries = 0, bool single_pulse = 0, int series = 13, int first = 0, int last = 10, const char* treename = "Waveforms", TString dir = "../Unpacked"){//"/projects/prospect/converted_data/Unpacked/DryCommissioning/"){

  time_t start, end;
  time(&start);
  
  Chain chain(series, first, last, dir, treename);
  chain.SetVerbose(1);
  chain.CreateChain();
  TChain *ch  = chain.GetChain();
  ULong64_t trig;//event/trigger number propagated from the input raw file, for
                 //tracking back to the original source,
  Long64_t noff;//time offset of first sample from start of file, in sample
                 //clock units (i.e.  4 ns steps forour digitizers)
  Int_t chan;//channel number (PMT number)
  UInt_t nsamp;//number of samples in the digitized waveform
  Short_t samps[kNmaxSamp];
  Int_t samples[kNmaxSamp], t[kNmaxSamp];
  Double_t samples_D[kNmaxSamp], t_D[kNmaxSamp], n[kNmaxSamp];
  memset(samples_D, 0, sizeof(samples_D[0])*kNmaxSamp);
  memset(n, 0, sizeof(n[0])*kNmaxSamp);
  ch->SetBranchAddress("trig", &trig);
  ch->SetBranchAddress("noff", &noff);
  ch->SetBranchAddress("chan", &chan);
  ch->SetBranchAddress("nsamp", &nsamp);
  ch->SetBranchAddress("samps", &samps);
  TGraph *gr;
  
  //create time/sample variable
  for(int i=0;i<kNmaxSamp;++i) {
    t[i] = i;
    t_D[i] = (double)i;
  }
  int skip = 0, nmax = 0;
  char sk[200];
  Long64_t nent = (max_entries == 0 ? ch->GetEntries() : max_entries);
  Long64_t mod = nent/10;
  cout<<ch->GetEntries()<<" entries\n";
  cout<<"|-----------------|"<<endl;
  for(int i=first; i<nent; ++i){
    if(i%mod==0){
      cout<<"* "<<flush;
      mod += nent/10;
    }
    ch->GetEntry(i);
    if((int)nsamp > nmax)nmax = (int)nsamp;
    
    if(single_pulse && chan != -100){
      for(ULong64_t j=0;j<nsamp;++j)samples[j] = (int)samps[j];
      gr = new TGraph(nsamp,t,samples);
      gr->SetTitle(Form("Channel %i",chan));
      gr->Draw("alp");
      gPad->Update();
      //sleep(1);
      cin>>sk;
      if(strcmp(sk,"q")==0)break;
    }else if(chan > -1 && chan < 308){
      for(ULong64_t j=0;j<nsamp;++j){
	++n[j];
	samples_D[j] = samples_D[j]/n[j]+ (n[j]-1)/n[j] * (double)samps[j];
      }
    }
  }
  cout<<endl;
  gr = new TGraph(nmax, t_D, samples_D);
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(0.6);
  gr->SetMarkerColor(kBlue);
  gr->SetLineColor(kRed);
  gr->Draw("alp");
  gr->SetTitle(Form("Average Pulse Waveform (Series %i)", series));
  gr->GetXaxis()->SetTitle("Sample Number");
  gr->GetYaxis()->SetTitle("ADC units");
  gPad->Update();
  time(&end);
  cout<<"Execution time: "<<difftime(end,start)<<" seconds for "<<nent<<" pulses"<<endl;
  return 0;
}
