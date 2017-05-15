//Start by creating a TChain of all useful runs. This can be done using
//CreateChain.C which includes all runs in the file runlist.dat. Make sure
//there are no runs from series018. It is also useful to cut out runs with
//high rates which are typically source calibration runs. One way of doing this
//is to compare the prompt and far windows to see what percentage of good events
//are accidental. This ratio is recorded in the last column of runlist.dat, so
//in CreateChain.C set the condition to add a run only if the ratio is low.
//Note that the far (time-displaced) window is 4 times longer than the prompt
//window.

{
  bool all = false;
  int plot = 1;
  TChain *ch = CreateChain();
  double mean = 0.88, rms = 0.056;
  TH1D *h;
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(1111);
  gStyle->SetTitleW(0.9);
  gStyle->SetStatX(0.695);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.25);
  TCanvas *c = new TCanvas("c","c",0,0,1000,600);
  //restrict to events in low evergy region with single prompt candidate 
  //and reconstructed inside cell volume   
  ch->Draw(">>list","prompt.E<4&&delayed.E<3&&prompt.E>0&&nPrompt==1&&abs(delayed.y)<600&&abs(prompt.y)<600");
  TEventList *list = (TEventList*)gDirectory->Get("list");
  ch->SetEventList(list);

  if(0){
    c->SetLogy();
    ch->Draw("delayed.E>>h(200)", "delayed.E<2");
    h = (TH1D*)gDirectory->Get("h");
    gPad->Update();
    h->SetTitle("Events near #alpha-candidate Energy Region (prompt energy > 0.2 MeV)");
    h->GetYaxis()->SetTitle("Number of Delayed Events");
    h->GetXaxis()->SetTitle("Delayed Energy (MeV)");
    TF1 *gs = new TF1("gs","gaus",0.82,0.98);
    gs->SetParameters(0.88,0.06);
    h->Fit(gs, "r");
    mean = gs->GetParameter(0);
    rms = gs->GetParameter(1);
    cout<<"mean: "<<mean<<"  rms: "<<rms<<endl;
    c->SaveAs("../plots/gaus_delayedE.png");
  }


  //Plot 1: 2D plot of prompt versus delayed energy
  if(plot==1 || all){
    ch->Draw("prompt.E:delayed.E>>h2(500,0,3,500,0,3,100)","nPrompt==1&&prompt.E<4&&prompt.E>0&&delayed.E<3","colz");
    TH2D *h2 = (TH2D*)gDirectory->Get("h2");
    gPad->Update();
    h2->SetTitle("Prompt versus Delayed Energy for BiPo Candidates");
    h2->GetYaxis()->SetTitle("Prompt Energy (MeV)");
    h2->GetXaxis()->SetTitle("Delayed Energy (MeV)");
    //    ch->Draw("prompt.E:delayed.E", Form("nPrompt==1&&prompt.E<4&&prompt.E>0.2&&delayed.E<0.94&&delayed.E>0.83","same");

  }
}

