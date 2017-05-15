{
TChain *ch = new TChain("BiPo");
gSystem->cd("/home/jonesdc/prospect/data/Analyzed/Phys_20161206/P50D/series007");
ch->Add("BiPoEvents_s007_f00000_ts1467830099_Phys.root");
ch->Add("BiPoEvents_s007_f00001_ts1467830551_Phys.root");
ch->Add("BiPoEvents_s007_f00002_ts1467830995_Phys.root");
ch->Add("BiPoEvents_s007_f00003_ts1467831428_Phys.root");
ch->Add("BiPoEvents_s007_f00004_ts1467831863_Phys.root");
ch->Add("BiPoEvents_s007_f00005_ts1467832291_Phys.root");
ch->Add("BiPoEvents_s007_f00006_ts1467832757_Phys.root");
ch->Add("BiPoEvents_s007_f00007_ts1467833256_Phys.root");
ch->Add("BiPoEvents_s007_f00008_ts1467833729_Phys.root");
ch->Add("BiPoEvents_s007_f00009_ts1467834278_Phys.root");
ch->Add("BiPoEvents_s007_f00010_ts1467834723_Phys.root");
ch->Add("BiPoEvents_s007_f00011_ts1467835181_Phys.root");
ch->Add("BiPoEvents_s007_f00012_ts1467835629_Phys.root");
ch->Add("BiPoEvents_s007_f00013_ts1467836072_Phys.root");
}
