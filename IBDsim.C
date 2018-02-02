{
  TChain *ch = ChainIBD(0,1199,"","IBD");
  TChain *chbg = ChainIBD(0,1199);
  TCanvas *c = new TCanvas("c","c",0,0,700,500);

  // ch->Draw("E>>h(100,0,9)","","goff");
  // TH1D *h = (TH1D*)gDirectory->Get("h");
  // h->SetLineColor(kBlue);
  // chbg->Draw("E>>hbg(100,0,9)","","goff");
  // TH1D *hbg = (TH1D*)gDirectory->Get("hbg");
  // hbg->SetLineColor(kRed);
  // hbg->Draw();
  // //  hbg->Scale((double)h->GetEntries()/(double)hbg->GetEntries());
  // h->Draw("same");
  c->SetLogy();
  //  //  TH1D *hd = (TH1D*)h->Clone();
  //  hd->Add(hbg,-1);
  //  hd->SetLineColor(kMagenta);
  //  hd->Draw("same");

}
