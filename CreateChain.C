#include <iostream>
#include <fstream>
#include "TChain.h"
#include "string.h"


TChain* CreateChain(){
  TChain *ch = new TChain("BiPo");
  //  ifstream f("goodruns.dat");
  ifstream f("runlist.dat");
  if(!f.is_open()){
    cout<<"Cannot open file. Exiting\n";
    return 0;
  }
  TString name;
  double a, b, c;
  string line;
  int n = 0;
  while(!f.eof()){
    f>>name>>a>>b>>c;
    getline(f, line);
    TString st(name(0,4));
    TString dir = "/home/jonesdc/prospect/data/Analyzed/Phys_20161206/P50D/serie"+ st + "/biPoEvents_" + name ;
    if(c>0.1)continue;
    ch->Add(dir.Data());
    cout<<dir.Data()<<endl;
    ++n;
    f.peek();
  }
  return ch;
}
