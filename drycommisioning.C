#include <iostream>
#include "TH2D.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"

int drycommissioning(int series = 10, double dt = 10 ){
  gStyle->SetOptStat(0);
  int nAdded = 0;
  double nSec = 10000;
  const int nFiles = 5;
  TChain *ch = new TChain("DetPulse");
  for(int i=0;i<nFiles;++i)
    nAdded += ch->Add(Form("/projects/prospect/converted_data/Pulse_Latest/DryCommissioning/series%03d/s%03d_f%05i*.root", series, series, i));

  printf("%i files added to chain. %lld total entries.\n", nAdded, ch->GetEntries());


  //  ch->Draw("det:t>>h(1000,-110,4.5e9,430,0,320","Entry$<100000&abs(5e7-t%5e7)<5e6");
  double pmt[14][11][5], N[14][11][5];
  for(int i=0;i<14;++i)
    for(int j=0;j<11;++j){
      for(int k=0;k<5;++k){
	pmt[i][j][k] = 0;
	N[i][j][k] = 0;
      }
    }

  Int_t det, col[5] = {kBlack, kRed, kBlue, kGreen+2, kCyan};
  TString title[5] = {"9.6 V", "10 V", "11 V", "12 V", "13 V"};
  Double_t t0, t_start=-1, t1 = 0, prev_time = -1, tx0, tx1;
  Int_t det;
  Long64_t evt;
  Float_t a, b, h, rise, PSD;
  Double_t t;
  ch->SetBranchAddress("evt", &evt);
  ch->SetBranchAddress("det", &det);
  ch->SetBranchAddress("t", &t);
  ch->SetBranchAddress("a", &a);
  ch->SetBranchAddress("b", &b);
  ch->SetBranchAddress("h", &h);
  ch->SetBranchAddress("rise", &rise);
  ch->SetBranchAddress("PSD", &PSD);
  TH1D *ha[nFiles];
  TH1D *hh[nFiles];
  TH1D *h1[nFiles];
  TH2D *h2 = new TH2D("h2", "h2", 1000, 0, 100, 350, -10, 340);
  TH1D *ht = new TH1D("ht","ht",100,-dt-5,dt+5);
  h2->SetCanExtend(TH1::kXaxis);
  printf("%lld entries.\n", ch->GetEntries());
  int nf = 0;
  h1[nf] =  new TH1D(Form("h1[%i]", nf), title[nf].Data(), 308,0,308);
  ha[nf] = new TH1D(Form("ha[%i]", nf), title[nf].Data(), 500, 0, 1000);
  hh[nf] = new TH1D(Form("hh[%i]", nf), title[nf].Data(), 500, 0, 500);
  for(int i=0;i<ch->GetEntries();++i){
    if(i%10000000==0)printf("Analyzing entry %i\n", i);
    ch->GetEntry(i);
    if(t<prev_time) {
      t1 += prev_time;
      ++nf;
      h1[nf] =  new TH1D(Form("h1[%i]", nf), title[nf].Data(), 308,0,308);
      ha[nf] = new TH1D(Form("ha[%i]", nf), title[nf].Data(), 500, 0, 1000);
      hh[nf] = new TH1D(Form("hh[%i]", nf), title[nf].Data(),500, 0, 500);
    }
    prev_time = t;
    t += t1;
    if(t/1e9 > nSec)break;
    if(det==-100){
      t0 = t;
      if(t_start==-1)t_start = t;
      //search backward for coincidences
      int n=0;
      double tpmt0 = 0, tpmt1 = 0;

      while(t0-t < dt){
	h1[nf]->Fill(det);
	h2->Fill(t/1e9, det);
	++n;
	if(i-n < 0) break;
	ch->GetEntry(i-n);
	t += t1;
	if(det==2 && tpmt0 == 0)tpmt0 = t;
	if(det==3 && tpmt1 == 0)tpmt1 = t;
	if(tpmt0 !=0 &&tpmt1 != 0)ht->Fill(tpmt0-tpmt1);
	if(n>310){
	  printf("Yikes\n");
	  break;
	}
	int x = int((det%28)/2), y = int(det/28);
	if(det>-1 && det<309){
	  pmt[x][y][nf] = pmt[x][y][nf]*N[x][y][nf]/(N[x][y][nf] + 1) + a/(N[x][y][nf] + 1);
	  ++N[x][y][nf];
	}
	ha[nf]->Fill(a);
	hh[nf]->Fill(h);
      }
      //search forward for coincidences
      ch->GetEntry(i);
      n=0;
      while(t-t0 < dt){
	++n;
	if(i+n>=ch->GetEntries()){
	  i+=n;
	  break;
	}
	ch->GetEntry(i+n);
	t += t1;
	if(det==2 && tpmt0 == 0)tpmt0 = t;
	if(det==3 && tpmt1 == 0)tpmt1 = t;
	if(tpmt0 !=0 && tpmt1 != 0)ht->Fill(tpmt0-tpmt1);
	if(n>310){
	  printf("Yikes\n");
	  break;
	}
	int x = int((det%28)/2), y = int(det/28);
	if(det>-1 && det<309){
	  pmt[x][y][nf] = pmt[x][y][nf]*N[x][y][nf]/(N[x][y][nf] + 1) + a/(N[x][y][nf] + 1);
	  ++N[x][y][nf];
	}
	h1[nf]->Fill(det);
	h2->Fill(t/1e9, det);
	ha[nf]->Fill(a);
	hh[nf]->Fill(h);
      }
      i += n-1;
    }
  }
  TCanvas *c = new TCanvas("c","c", 0, 0, 1800,1000);
  c->Divide(2,2);
  c->cd(1);

  h2->SetTitle(Form("Signals  within %0.0f ns of OCS signal (Series %i)", dt, series));
  h2->GetYaxis()->SetTitle("PMT #");
  h2->GetXaxis()->SetTitle("Time (s)");
  h2->Draw("colz");
  TPaveText *tt[5];
  for(int i=0;i<=nf;++i){
    tt[i] = new TPaveText(0.102 + i*0.15, 0.86,0.1+0.15*(i+1),0.898,"ndc");
    tt[i]->AddText(title[i].Data());
    tt[i]->SetFillColor(0);
    tt[i]->SetBorderSize(0);
    tt[i]->SetShadowColor(0);
    tt[i]->SetTextColor(kRed);
    tt[i]->Draw();
  }
  c->cd(2);
  ha[0]->SetTitle("Pulse Integral Spectrum of Signals Coincident with OCS");
  ha[0]->GetYaxis()->SetLimits(0,  ha[0]->GetMaximum()*1.1);
  ha[0]->GetYaxis()->SetRangeUser(0,  ha[0]->GetMaximum()*1.1);
  ha[0]->GetXaxis()->SetTitle("Pulse Integral (ADC units)");
  ha[0]->SetLineColor(col[0]);
  ha[0]->SetLineWidth(2);
  ha[0]->Draw();
  cout<<nf<<endl;
  TPaveText *pt = new TPaveText(0.85,0.55,0.99,0.9,"ndc");
  pt->SetFillColor(0);
  pt->SetShadowColor(0);
  pt->AddText(Form("- %s", title[0].Data()));
  ((TText*)pt->GetListOfLines()->Last())->SetTextColor(col[0]);
  for(int i=1;i<=nf;++i){
    pt->AddText(Form("- %s", title[i].Data()));
    ((TText*)pt->GetListOfLines()->Last())->SetTextColor(col[i]);
    ha[i]->SetLineColor(col[i]);
    ha[i]->SetLineWidth(2);
    ha[i]->Draw("same");
  }
  pt->Draw();
  c->cd(3);
  h1[0]->SetTitle(Form("Signals within %0.0f ns of OCS signal (Series %i)", dt, series));
  h1[0]->GetXaxis()->SetTitle("PMT #");
  h1[0]->SetLineColor(col[0]);
  h1[0]->Draw();
  for(int i=1;i<=nf;++i){
    h1[i]->SetLineColor(col[i]);
    h1[i]->Draw("same");
  }
  pt->Draw();

  c->cd(4);
  hh[0]->SetTitle("Pulse Height Spectrum of Signals Coincident with OCS");
  hh[0]->GetXaxis()->SetTitle("Pulse Integral (ADC units)");
  hh[0]->SetLineColor(col[0]);
  hh[0]->SetLineWidth(2);
  hh[0]->Draw();

  for(int i=1;i<=nf;++i){
    hh[i]->SetLineColor(col[i]);
    hh[i]->SetLineWidth(2);
    hh[i]->Draw("same");
  }
  pt->Draw();
  printf("Total time %f seconds.\n", (t-t_start)*1e-9);
  c->SaveAs("../plots/dry_commissioning_plots.pdf");
  c->SaveAs("../plots/dry_commissioning_plots.png");


  gStyle->SetPadRightMargin(0.2);
  gStyle->SetTitleW(0.9);
  TCanvas *c1 = new TCanvas("c1","c1", 0, 0, 1500,1000);
  c1->Divide(3,2);
  TCanvas *c2 = new TCanvas("c2","c2", 0, 0, 1500,1000);
  c2->Divide(3,2);
  TCanvas *c3 = new TCanvas("c3","c3", 0, 0, 1500,1000);
  c3->Divide(3,2);
  TH2D *hN[5];
  TH2D *hdet[5];
  for(int nfi=0;nfi<=nf;++nfi){
    hN[nfi] = new TH2D(Form("hN[%i]",nfi), Form("hN[%i]",nfi),14,0,14,11,0,11);
    hN[nfi]->SetTitle(Form("Number of Pulses Coincident with OCS by Cell, (Series %i, V=%s)",series, title[nfi].Data()));
    hN[nfi]->GetXaxis()->SetTitle("Column");
    hN[nfi]->GetYaxis()->SetTitle("Row");
    hdet[nfi] = new TH2D(Form("hdet[%i]",nfi),Form("hdet[%i]",nfi),14,0,14,11,0,11);
    hdet[nfi]->SetTitle(Form("Average Pulse Integral by Cell (Series %i, V=%s)", series, title[nfi].Data()));
    hdet[nfi]->GetXaxis()->SetTitle("Column");
    hdet[nfi]->GetYaxis()->SetTitle("Row");
    for(int i=0;i<14;++i){
      for(int j=0;j<11;++j){
	hN[nfi]->Fill(i,j,N[i][j][nfi]);
	hdet[nfi]->Fill(i,j,pmt[i][j][nfi]);
      }
    }
    c1->cd(nfi+1);
    hN[nfi]->Draw("colz");
    c2->cd(nfi+1);
    hdet[nfi]->Draw("colz");
    c3->cd(nfi+1)->SetLogz();
    hdet[nfi]->Draw("colz");
  }
  c2->SaveAs("../plots/dry_commissioning_OCS_cell_average_integral.png");
  c1->SaveAs("../plots/dry_commissioning_OCS_npulses.png");
  c3->SaveAs("../plots/dry_commissioning_OCS_cell_average_integral_logz.png");
  return 0;
}

