// builds runlist from rtree.root

void make_run_list_consistent(const char * infile="rtree.root")
{
  TFile * tf = new TFile(infile,"READ");
  TTree * tr = (TTree*) tf->Get("rellum");
  tr->SetScanField(10000);
  char outfile[64];
  strcpy(outfile,"run_table_consistent.txt");
  gSystem->RedirectOutput(outfile,"w");
  tr->Scan("i:runnum:fi:fill","isConsistent");
  gSystem->RedirectOutput(0);
  printf("%s created\n",outfile);
};
