// shows some plots based on different scaler bit combinations
// - total counts vs. scaler bit (bar chart)
// - total counts / bXings vs. bXing
//
// (see rellumi.pdf)
//

void bit_combinations(const char * filename="counts.root")
{
  TFile * infile = new TFile(filename,"READ");
  TTree * acc = (TTree*) infile->Get("acc");

  // set branch addresses
  Int_t index,runnum,fill_index,fill;
  Double_t t;
  Double_t bbc[8];
  Double_t zdc[8];
  Double_t vpd[4];
  Double_t tot_bx;
  Int_t blue,yell;
  acc->SetBranchAddress("i",&index);
  acc->SetBranchAddress("runnum",&runnum);
  acc->SetBranchAddress("fi",&fill_index);
  acc->SetBranchAddress("fill",&fill);
  acc->SetBranchAddress("t",&t);
  acc->SetBranchAddress("tot_bx",&tot_bx);
  acc->SetBranchAddress("blue",&blue);
  acc->SetBranchAddress("yell",&yell);
  char bbc_br[8][16];
  char zdc_br[8][16];
  char vpd_br[4][16];
  for(Int_t i=0; i<8; i++)
  {
    sprintf(bbc_br[i],"bbc_%d",i);
    sprintf(zdc_br[i],"zdc_%d",i);
    acc->SetBranchAddress(bbc_br[i],&(bbc[i]));
    acc->SetBranchAddress(zdc_br[i],&(zdc[i]));
    if(i<4) 
    {
      sprintf(vpd_br[i],"vpd_%d",i);
      acc->SetBranchAddress(vpd_br[i],&(vpd[i]));
    };
  };

  // ----------------------------------------------

  // total counts vs. scaler bit bar charts
  TH1F * ntot_vs_bits_bbc = new TH1F()
  TH1F * ntot_vs_bits_zdc = new TH1F()
  TH1F * ntot_vs_bits_vpd = new TH1F()


  // scaler bit combination names
  char comb[8][16];
  strcpy(comb[0],"none");
  strcpy(comb[1],"E");
  strcpy(comb[2],"W");
  strcpy(comb[3],"W+E");
  strcpy(comb[4],"X");
  strcpy(comb[5],"X+E");
  strcpy(comb[6],"X+W");
  strcpy(comb[7],"X+W+E");

  // fill bar charts
  for(Int_t i=0; i<acc->GetEntries(); i++)
  {
    acc->GetEntry(i);
    for(Int_t j=0; j<8; j++)
    {
      ntot_vs_bits_bbc->Fill(comb[j],bbc[j]);
      ntot_vs_bits_zdc->Fill(comb[j],zdc[j]);
      if(j<4) ntot_vs_bits_vpd->Fill(comb[j],vpd[j]);
    };
  };

  ntot_vs_bits_bbc->SetStats(0);
  ntot_vs_bits_bbc->SetTitle("total bbc counts vs. scaler bits");
  ntot_vs_bits_bbc->SetBarWidth(0.4);
  ntot_vs_bits_bbc->SetBarOffset(0.55);
  ntot_vs_bits_bbc->SetFillColor(50);
  TCanvas * c_bbc_bits = new TCanvas("c_bbc_bits","c_bbc_bits",700,500);
  ntot_vs_bits_bbc->Draw("bar2");

  ntot_vs_bits_zdc->SetStats(0);
  ntot_vs_bits_zdc->SetTitle("total zdc counts vs. scaler bits");
  ntot_vs_bits_zdc->SetBarWidth(0.4);
  ntot_vs_bits_zdc->SetBarOffset(0.55);
  ntot_vs_bits_zdc->SetFillColor(50);
  TCanvas * c_zdc_bits = new TCanvas("c_zdc_bits","c_zdc_bits",700,500);
  ntot_vs_bits_zdc->Draw("bar2");

  ntot_vs_bits_vpd->SetStats(0);
  ntot_vs_bits_vpd->SetTitle("total vpd counts vs. scaler bits");
  ntot_vs_bits_vpd->SetBarWidth(0.4);
  ntot_vs_bits_vpd->SetBarOffset(0.55);
  ntot_vs_bits_vpd->SetFillColor(50);
  TCanvas * c_vpd_bits = new TCanvas("c_vpd_bits","c_vpd_bits",700,500);
  ntot_vs_bits_vpd->Draw("bar2");

  // ----------------------------------------------


  // get maximum number of runs  = IMAX
  Int_t IMAX_tmp = acc->GetMaximum("i");
  const Int_t IMAX = IMAX_tmp;


  // compute no. bXings with possible interaction in a given run
  Int_t nbx[IMAX]; 
  for(Int_t i=0; i<IMAX; i++) nbx[i]=0;
  for(Int_t i=0; i<acc->GetEntries(); i++)
  {
    acc->GetEntry(i);
    if(blue*yell != 0) nbx[index-1]++;
  };
};





