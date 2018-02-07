
#include <ctime>
#include "Chain.cc"
#include "TF1.h"
#include "TFitResult.h"
#include "TGaxis.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <vector>
#include <string.h>

using namespace std;
const int kNbinsX = 200;
const double kMaxDT = 1000; //maximum time (ns) between coincident pulses
const double kMaxSPE = 400; //maximum ADC value included in SPE histograms

int CrossTalk(int npulses = 0, int series = 13, int first = 0, int last = 0, const char* treename = "DetPulse", TString dir = "../DetPulse"){//TString dir = "/projects/prospect/converted_data/Pulse_Latest/DryCommissioning/"){
  time_t start, end;
  time(&start);
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.9);
  gStyle->SetTitleW(0.9);
  Chain chain(series, first, last, dir, treename);
  chain.SetVerbose(1);
  if(npulses>0)
    chain.SetMaxPulses(npulses);
  chain.CreateChain();
  TChain *ch  = chain.GetChain();
  Int_t det;
  Long64_t evt, prev_evt = 0;
  Float_t a, b, h, rise, PSD;
  Double_t t;
  ch->SetBranchStatus("*",0);
  //ch->SetBranchStatus("evt", 1);
  ch->SetBranchStatus("det", 1);
  ch->SetBranchStatus("a", 1);
  ch->SetBranchStatus("t", 1);
  //ch->SetBranchAddress("evt", &evt);
  ch->SetBranchAddress("det", &det);
  ch->SetBranchAddress("t", &t);
  ch->SetBranchAddress("a", &a);
  // ch->SetBranchAddress("b", &b);
  // ch->SetBranchAddress("h", &h);
  // ch->SetBranchAddress("rise", &rise);
  // ch->SetBranchAddress("PSD", &PSD);
  int nx = 14, ny = 11, nc = nx*ny;
  int *celle = new int[nc],*cellw = new int[nc];
  double first_time, prev_time = 0;
  memset(celle, 0, sizeof(celle[0])*nc);
  memset(cellw, 0, sizeof(cellw[0])*nc);
  double *te = new double[nc], *tw = new double[nc];
  memset(te, 0, sizeof(te[0])*nc);
  memset(tw, 0, sizeof(tw[0])*nc);
  vector<int>vCell;
  vector<double>vArea;
  vector<double>vDt;
  TH1D *he[16], *hw[16], *hSPEe[16], *hSPEw[16];
  TH1D *hdt = new TH1D("hdt","Time of Pulses Relative to OCS Pulse",kMaxDT,-2*kMaxDT,2*kMaxDT);
  hdt->SetXTitle("#Delta t (ns)");
  TH1D *hdtn = new TH1D("hdtn","hdtn",kMaxDT,-2*kMaxDT,2*kMaxDT);
  hdtn->SetLineColor(kRed);
  TH1D *hdtp = new TH1D("hdtp","hdtp",kMaxDT,-2*kMaxDT,2*kMaxDT);
  hdtp->SetLineColor(kGreen+2);
  Double_t binEdge[kNbinsX+1];
  Double_t lBin = 10, hBin = 100000, log10l = log10(lBin), log10h = log10(hBin);
  for(int i=0;i<=kNbinsX;++i){
    binEdge[i] = pow(10,log10l + i*(log10h-log10l)/kNbinsX);
    //cout<<"binEdge["<<i<<"]: "<<binEdge[i]<<endl;
  }
  for(int i=0;i<16;++i){
    hSPEe[i] = new TH1D(Form("hSPEe%i",i),Form("hSPEe%i",i),200,0,kMaxSPE);
    hSPEe[i]->SetLineColor(kRed);
    hSPEw[i] = new TH1D(Form("hSPEw%i",i),Form("hSPEw%i",i),200,0,kMaxSPE);
    hSPEw[i]->SetLineColor(kBlue);
    he[i] = new TH1D(Form("he[%i]",i),Form("PMT Integral, Cell X=%i, Y=%i",i%4+5,i/4+3),kNbinsX, binEdge);
    he[i]->SetLineColor(kRed);
    he[i]->GetXaxis()->SetTitle("ADC Channels");
    he[i]->GetXaxis()->SetTitleOffset(1.3);
    he[i]->GetXaxis()->SetTitleSize(0.6);
    hw[i] = new TH1D(Form("hw[%i]",i),Form("Cell X=%i, Y=%i Integral",i%4+5,i/4+3),kNbinsX,binEdge);
    hw[i]->SetLineColor(kBlue);
    hw[i]->GetXaxis()->SetTitle("ADC Channels");
    hw[i]->GetXaxis()->SetTitleOffset(1.3);
    hw[i]->GetXaxis()->SetTitleSize(0.6);

  }
  TH2D *h2e = new TH2D("h2e","PMT Hit Map East",nx,0,nx,ny,0,ny);
  h2e->SetXTitle("Column");
  h2e->SetYTitle("Row");
  TH2D *h2w = new TH2D("h2w","PMT Hit Map West",nx,0,nx,ny,0,ny);
  h2w->SetXTitle("Column");
  h2w->SetYTitle("Row");
  TH2D *h2c = new TH2D("h2c","PMT Hit Map Cell Coincidence",nx,0,nx,ny,0,ny);
  TH2D *h2dt = new TH2D("h2dt","Cell #Delta t",nc, 0, nc, 200, -2*kMaxDT, +2*kMaxDT);
  //  TH2D *h2dt = new TH2D("h2dt","h2dt",nx,0,nx,ny,0,ny);
  h2c->SetXTitle("Column");
  h2c->SetYTitle("Row");
  int coinc = 0;
  Long64_t nent = (npulses == 0 ? ch->GetEntries() : npulses);
  Long64_t mod = nent/10;
  cout<<"|-----------------|"<<endl;
  ch->GetEntry(0);
  first_time = t;
  for(int i=0;i<nent;++i){
    if(i%mod==0){
      cout<<"* "<<flush;
      mod += nent/10;
    }
    ch->GetEntry(i);
    if(det == -100){
      int ent = 0;
      //Look backward for coincidences
      double t_ocs = t;
      while(i-ent > -1){
	ch->GetEntry(i-ent);
	double dt = t_ocs - t;
	int cl = det/2;
	int x = int(det/2)%nx, y = int(det/2)/nx, idx = x-5+(y-3)*4;
	if(dt > kMaxDT) break;
	if(det<0 || det > 2*nx*ny -1){
	  ++ent;
	  continue;
	}
	   
	//Fill histograms here
	if(det%2==0){
	  h2e->Fill(x,y);
	  if(x>4 && x<9 && y>2 && y<7){
	    hdt->Fill(-dt);
	    if(x>5&&x<8&&y>3&&y<6)
	      hdtp->Fill(-dt);
	    else
	      hdtn->Fill(-dt);
	    he[idx]->Fill(a);
	    if(a < kMaxSPE)
	      hSPEe[idx]->Fill(a);
	  }
	}else{
	  h2w->Fill(x,y);
	  if(x>4 && x<9 && y>2 && y<7){
	    hdt->Fill(-dt);
	    if(x>5&&x<8&&y>3&&y<6)
	      hdtp->Fill(-dt);
	    else
	      hdtn->Fill(-dt);
	    hw[idx]->Fill(a);
	    if(a < kMaxSPE)
	      hSPEw[idx]->Fill(a);
	  }
	}
	if(det%2==0){
	  ++celle[cl];
	  te[cl] = t;
	}else{
	  ++cellw[cl];
	  tw[cl] = t;
	}
	if(cellw[cl]>0&&celle[cl]==cellw[cl]){
	  vCell.push_back(cl);
	  vDt.push_back(te[cl]-tw[cl]);
	}
	++ent;
      }
      
      //Look forward for coincidences
      ent = 1;
      while(i+ent < ch->GetEntries()){
	if((i+ent)%mod==0){
	  cout<<"* "<<flush;
	  mod += nent/10;
	}

	ch->GetEntry(i+ent);
	double dt = t - t_ocs;
	int cl = det/2;
	int x = int(det/2)%nx, y = int(det/2)/nx, idx = x-5+(y-3)*4;
	if(dt > kMaxDT) break;
	if(det<0 || det > 2*nx*ny -1){
	  ++ent;
	  continue;
	}
	
	
	//Fill histograms here
	if(det%2==0){
	  h2e->Fill(x,y);
	  if(x>4 && x<9 && y>2 && y<7){
	    hdt->Fill(dt);
	    if(x>5&&x<8&&y>3&&y<6)
	      hdtp->Fill(dt);
	    else
	      hdtn->Fill(dt);
	    he[idx]->Fill(a);
	    if(a < kMaxSPE)
	      hSPEe[idx]->Fill(a);
	  }
	}else{
	  h2w->Fill(x,y);
	  if(x>4 && x<9 && y>2 && y<7){
	    hdt->Fill(dt);
	    if(x>5&&x<8&&y>3&&y<6)
	      hdtp->Fill(dt);
	    else
	      hdtn->Fill(dt);
	    hw[idx]->Fill(a);
	    if(a < kMaxSPE)
	      hSPEw[idx]->Fill(a);
	  }
	}
	if(det%2==0){
	  ++celle[cl];
	  te[cl] = t;
	}else{
	  ++cellw[cl];
	  tw[cl] = t;
	}
	if(cellw[cl]>0&&celle[cl] == cellw[cl]){
	  vCell.push_back(cl);
	  vDt.push_back(te[cl] - tw[cl]);
	}
	++ent;
      }
      i += ent - 1;
    }

    for(int c=0;c<int(vCell.size());++c){
      int x = vCell[c]%nx, y = vCell[c]/nx;
      h2c->Fill(x,y);
      h2dt->Fill(vCell[c], vDt[c]);
    }
    vCell.resize(0);
    vDt.resize(0);
    memset(te, 0, sizeof(te[0])*nc);
    memset(tw, 0, sizeof(tw[0])*nc);
    memset(celle, 0, sizeof(celle[0])*nc);
    memset(cellw, 0, sizeof(cellw[0])*nc);
    prev_evt = evt;
    prev_time = t;
  }

  printf("\n");

  TF1 *f = new TF1("f","[0]*exp(-pow((x-[1])/(2*sqrt([2]*[2]+[4]*[4])),2))+[3]*exp(-pow((x-2*[1])/(2*sqrt(2*[2]*[2]+[4]*[4])),2))",0,250);
  double spee[16], spew[16];
  TGaxis *axise[16], *axisw[16];
  double lLim = 30,hLim = 80;
  for(int i=0;i<16;++i){
    double height1 = hSPEe[i]->GetMaximum(), height2 = height1/10.0;
    double mean = hSPEe[i]->GetBinCenter(hSPEe[i]->GetMaximumBin());
    if(mean>hLim){
      mean/=2;
      height1 = height2;
      height2 *= 10;
    }
    f->SetParameters(height1, mean, sqrt(mean), height2, 0);
    // f->SetParLimits(0,0,1.1*hSPEe[i]->GetMaximum());
    // f->SetParLimits(1,lLim,hLim);
    // f->SetParLimits(2,0,1e6);
    // f->SetParLimits(3,0,hSPEe[i]->GetMaximum());
    // f->SetParLimits(4,0,1e6);
    // f->SetParLimits(5,0,hSPEe[i]->GetMaximum());
    f->SetLineColor(kRed);
    //f->SetParLimits(6,0,1e6);
    cout<<i<<"-------\n";
    if(!(hSPEe[i]->Fit("f","rs").Get()->IsValid()))hSPEe[i]->Fit("f","r");
    cout<<f->GetParameter(1)<<endl;
    spee[i] = f->GetParameter(1);
    height1 = hSPEw[i]->GetMaximum(), height2 = height1/10.0;
    mean = hSPEw[i]->GetBinCenter(hSPEw[i]->GetMaximumBin());    
    if(mean>hLim){
      mean/=2;
      height1 = height2;
      height2 *= 10;
    }
    f->SetParameters(height1, mean, sqrt(mean), height2, 0);
    // f->SetParLimits(0,0,1.1*hSPEw[i]->GetMaximum());
    // f->SetParLimits(3,0,hSPEw[i]->GetMaximum());
    // f->SetParLimits(5,0,hSPEw[i]->GetMaximum());
    f->SetLineColor(kBlue);
    if(!(hSPEw[i]->Fit("f","rs").Get()->IsValid()))hSPEw[i]->Fit("f","r");
    cout<<f->GetParameter(1)<<endl;
    cout<<"*******\n";

    spew[i] = f->GetParameter(1);
    axisw[i] = new TGaxis(0.1, 0.9, 0.99, 0.9, hw[i]->GetXaxis()->GetXmin()/spew[i],hw[i]->GetXaxis()->GetXmax()/spew[i], 510,"-");
  }

  TCanvas *c = new TCanvas("c","c",0,0,1600,1000);
  c->Divide(2,2);
  c->cd(1)->SetLogz();
  h2w->Draw("colz");
  c->cd(2)->SetLogz();
  h2e->Draw("colz");
  c->cd(3)->SetLogz();
  h2c->Draw("colz");
  c->cd(4);
  h2dt->Draw("colz");
  c->ForceUpdate();
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetLabelSize(0.06);
  gStyle->SetLabelSize(0.06,"Y");
  gStyle->SetPadTopMargin(0.08);

  int can[16] = {13,14,15,16,9,10,11,12,5,6,7,8,1,2,3,4};
  TCanvas *cSPEe = new TCanvas("cSPEe","cSPEe",0,0,1600,1000);
  cSPEe->Divide(4,4);
  TCanvas *cSPEw = new TCanvas("cSPEw","cSPEw",0,0,1600,1000);
  cSPEw->Divide(4,4);
  for(int i=0;i<16;++i){
    cSPEe->cd(can[i]);
    hSPEe[i]->Draw();
    cSPEw->cd(can[i]);
    hSPEw[i]->Draw();
  }
  TCanvas *c2 = new TCanvas("c2","c2",0,0,1600,1000);
  c2->Divide(4,4);
  double primary = 0, secondary = 0;
  for(int i=0;i<16;++i){
    c2->cd(can[i])->SetLogy();
    gPad->SetLogx();
    he[i]->Draw();
    gPad->Update();
    TAxis *axX =  he[i]->GetXaxis();
    double xmax = axX->GetXmax();
    double xmin = axX->GetXmin();
    double y =  4*he[i]->GetMaximum();
    he[i]->SetMaximum(y);
    y *= 0.6;
    //cout<<xmin<<" "<<xmax<<" "<<y<<" "<<xmin/spee[i]<<" "<<xmax/spee[i]<<endl;
    axise[i] = new TGaxis(xmin,y,xmax,y, xmin/spee[i],xmax/spee[i], 510,"GB");
    axise[i]->SetLabelSize(0.06);
    axise[i]->SetLineColor(kRed);
    //axise[i]->SetLabelColor(kRed);
    axise[i]->Draw();
    hw[i]->Draw("sames");
    //cout<<xmin<<" "<<xmax<<" "<<y<<" "<<xmin/spew[i]<<" "<<xmax/spew[i]<<endl;
    axisw[i] = new TGaxis(xmin,y,xmax,y, xmin/spew[i],xmax/spew[i], 510,"-GU");
    axisw[i]->SetLineColor(kBlue);
    //axisw[i]->SetLabelOffset(0.01);
    //    axisw[i]->SetLabelColor(kRed);
    axisw[i]->Draw();
    double integral = he[i]->Integral("width")/spee[i]+hw[i]->Integral("width")/spew[i];
    if(i==5||i==6||i==9||i==10){
      primary += integral;
    }else {
      secondary += integral;
    }
    cout<<"Integral: "<<integral<<endl;
  }
  cout<<"Primary: "<<primary<<endl;
  cout<<"Secondary: "<<secondary<<endl;
  cout<<secondary/primary*100<<"% signal leakage\n";
  TCanvas *c3 = new TCanvas("c3","c3",0,0,800,550);
  c3->SetLogy();
  hdt->Draw();
  hdtp->Draw("same");
  hdtn->Draw("same");
  hdt->GetXaxis()->SetTitleSize(0.05);
  gPad->Update();
  //  hSPEe[6]->Draw();
  //hSPEw[6]->Draw("sames");
  // double integral = 0;
  // for(int i=1;i<=he[0]->GetNbinsX();++i){
  //   integral += he[0]->GetBinContent(i)*he[0]->GetBinWidth(i);
  // }
  // cout<<"Int: "<<integral<<" "<<he[0]->Integral()<<" "<<he[0]->Integral("width")<<endl;
  for(int i=0;i<16;++i){
    cout<<"SPE east "<<i<<": "<<spee[i]<<endl;
    cout<<"SPE west "<<i<<": "<<spew[i]<<endl;
  }
  time(&end);
  cout<<"Execution time: "<<difftime(end,start)<<" seconds for "<<nent<<" pulses"<<endl;
  return 0;
}
