// combines sums.root and rdat_i.root into one final tree, with relative luminosities
// and total scaler counts for each run
//
// output is "rtree.root"

void combineAll(const char * rdatfile="rdat_i.root", const char * sumfile="sums.root")
{
  // open trees
  TFile * rdatfile_tf = new TFile(rdatfile,"READ");
  TTree * rdat = (TTree*) rdatfile_tf->Get("rtr");
  TFile * sumfile_tf = new TFile(sumfile,"READ");
  TTree * sums = (TTree*) sumfile_tf->Get("sum");


  // check if trees have same size
  if(rdat->GetEntries() != sums->GetEntries())
  {
    fprintf(stderr,"ERROR: tree sizes unequal\n");
    return;
  };

  Int_t i_rdat,runnum_rdat,fill_rdat;
  Double_t t_rdat;
  Float_t RRe[3][10]; // [tbit (bbc,zdc,vpd)] [rellum]
  Float_t RRw[3][10]; 
  Float_t RRx[3][10]; 
  Float_t RR[3][10]; // mean R
  Float_t RR_err[3][10]; // error of mean R
  Int_t i_sums,runnum_sums,fill_sums,fi_sums; // (fi not in rdat tree)
  Double_t t_sums,tau;
  Int_t num_runs;
  Double_t bbce,bbcw,bbcx;
  Double_t zdce,zdcw,zdcx;
  Double_t vpde,vpdw,vpdx;
  Double_t tot_bx;
  Int_t pattern;
  Float_t d_vz[3]; // [cbit]
  Float_t d_xx[3][3]; // [xbit] [tbit]
  Bool_t isConsistent; // true for a run if diagnostics are "consistent"

  rdat->SetBranchAddress("i",&i_rdat);
  rdat->SetBranchAddress("runnum",&runnum_rdat);
  rdat->SetBranchAddress("fill",&fill_rdat);
  rdat->SetBranchAddress("t",&t_rdat);

  rdat->SetBranchAddress("R1_bbce",&(RRe[0][1]));
  rdat->SetBranchAddress("R2_bbce",&(RRe[0][2]));
  rdat->SetBranchAddress("R3_bbce",&(RRe[0][3]));
  rdat->SetBranchAddress("R4_bbce",&(RRe[0][4]));
  rdat->SetBranchAddress("R5_bbce",&(RRe[0][5]));
  rdat->SetBranchAddress("R6_bbce",&(RRe[0][6]));
  rdat->SetBranchAddress("R7_bbce",&(RRe[0][7]));
  rdat->SetBranchAddress("R8_bbce",&(RRe[0][8]));
  rdat->SetBranchAddress("R9_bbce",&(RRe[0][9]));
  rdat->SetBranchAddress("R1_zdce",&(RRe[1][1]));
  rdat->SetBranchAddress("R2_zdce",&(RRe[1][2]));
  rdat->SetBranchAddress("R3_zdce",&(RRe[1][3]));
  rdat->SetBranchAddress("R4_zdce",&(RRe[1][4]));
  rdat->SetBranchAddress("R5_zdce",&(RRe[1][5]));
  rdat->SetBranchAddress("R6_zdce",&(RRe[1][6]));
  rdat->SetBranchAddress("R7_zdce",&(RRe[1][7]));
  rdat->SetBranchAddress("R8_zdce",&(RRe[1][8]));
  rdat->SetBranchAddress("R9_zdce",&(RRe[1][9]));
  rdat->SetBranchAddress("R1_vpde",&(RRe[2][1]));
  rdat->SetBranchAddress("R2_vpde",&(RRe[2][2]));
  rdat->SetBranchAddress("R3_vpde",&(RRe[2][3]));
  rdat->SetBranchAddress("R4_vpde",&(RRe[2][4]));
  rdat->SetBranchAddress("R5_vpde",&(RRe[2][5]));
  rdat->SetBranchAddress("R6_vpde",&(RRe[2][6]));
  rdat->SetBranchAddress("R7_vpde",&(RRe[2][7]));
  rdat->SetBranchAddress("R8_vpde",&(RRe[2][8]));
  rdat->SetBranchAddress("R9_vpde",&(RRe[2][9]));

  rdat->SetBranchAddress("R1_bbcw",&(RRw[0][1]));
  rdat->SetBranchAddress("R2_bbcw",&(RRw[0][2]));
  rdat->SetBranchAddress("R3_bbcw",&(RRw[0][3]));
  rdat->SetBranchAddress("R4_bbcw",&(RRw[0][4]));
  rdat->SetBranchAddress("R5_bbcw",&(RRw[0][5]));
  rdat->SetBranchAddress("R6_bbcw",&(RRw[0][6]));
  rdat->SetBranchAddress("R7_bbcw",&(RRw[0][7]));
  rdat->SetBranchAddress("R8_bbcw",&(RRw[0][8]));
  rdat->SetBranchAddress("R9_bbcw",&(RRw[0][9]));
  rdat->SetBranchAddress("R1_zdcw",&(RRw[1][1]));
  rdat->SetBranchAddress("R2_zdcw",&(RRw[1][2]));
  rdat->SetBranchAddress("R3_zdcw",&(RRw[1][3]));
  rdat->SetBranchAddress("R4_zdcw",&(RRw[1][4]));
  rdat->SetBranchAddress("R5_zdcw",&(RRw[1][5]));
  rdat->SetBranchAddress("R6_zdcw",&(RRw[1][6]));
  rdat->SetBranchAddress("R7_zdcw",&(RRw[1][7]));
  rdat->SetBranchAddress("R8_zdcw",&(RRw[1][8]));
  rdat->SetBranchAddress("R9_zdcw",&(RRw[1][9]));
  rdat->SetBranchAddress("R1_vpdw",&(RRw[2][1]));
  rdat->SetBranchAddress("R2_vpdw",&(RRw[2][2]));
  rdat->SetBranchAddress("R3_vpdw",&(RRw[2][3]));
  rdat->SetBranchAddress("R4_vpdw",&(RRw[2][4]));
  rdat->SetBranchAddress("R5_vpdw",&(RRw[2][5]));
  rdat->SetBranchAddress("R6_vpdw",&(RRw[2][6]));
  rdat->SetBranchAddress("R7_vpdw",&(RRw[2][7]));
  rdat->SetBranchAddress("R8_vpdw",&(RRw[2][8]));
  rdat->SetBranchAddress("R9_vpdw",&(RRw[2][9]));

  rdat->SetBranchAddress("R1_bbcx",&(RRx[0][1]));
  rdat->SetBranchAddress("R2_bbcx",&(RRx[0][2]));
  rdat->SetBranchAddress("R3_bbcx",&(RRx[0][3]));
  rdat->SetBranchAddress("R4_bbcx",&(RRx[0][4]));
  rdat->SetBranchAddress("R5_bbcx",&(RRx[0][5]));
  rdat->SetBranchAddress("R6_bbcx",&(RRx[0][6]));
  rdat->SetBranchAddress("R7_bbcx",&(RRx[0][7]));
  rdat->SetBranchAddress("R8_bbcx",&(RRx[0][8]));
  rdat->SetBranchAddress("R9_bbcx",&(RRx[0][9]));
  rdat->SetBranchAddress("R1_zdcx",&(RRx[1][1]));
  rdat->SetBranchAddress("R2_zdcx",&(RRx[1][2]));
  rdat->SetBranchAddress("R3_zdcx",&(RRx[1][3]));
  rdat->SetBranchAddress("R4_zdcx",&(RRx[1][4]));
  rdat->SetBranchAddress("R5_zdcx",&(RRx[1][5]));
  rdat->SetBranchAddress("R6_zdcx",&(RRx[1][6]));
  rdat->SetBranchAddress("R7_zdcx",&(RRx[1][7]));
  rdat->SetBranchAddress("R8_zdcx",&(RRx[1][8]));
  rdat->SetBranchAddress("R9_zdcx",&(RRx[1][9]));
  rdat->SetBranchAddress("R1_vpdx",&(RRx[2][1]));
  rdat->SetBranchAddress("R2_vpdx",&(RRx[2][2]));
  rdat->SetBranchAddress("R3_vpdx",&(RRx[2][3]));
  rdat->SetBranchAddress("R4_vpdx",&(RRx[2][4]));
  rdat->SetBranchAddress("R5_vpdx",&(RRx[2][5]));
  rdat->SetBranchAddress("R6_vpdx",&(RRx[2][6]));
  rdat->SetBranchAddress("R7_vpdx",&(RRx[2][7]));
  rdat->SetBranchAddress("R8_vpdx",&(RRx[2][8]));
  rdat->SetBranchAddress("R9_vpdx",&(RRx[2][9]));

  rdat->SetBranchAddress("R1_bbc_mean",&(RR[0][1]));
  rdat->SetBranchAddress("R2_bbc_mean",&(RR[0][2]));
  rdat->SetBranchAddress("R3_bbc_mean",&(RR[0][3]));
  rdat->SetBranchAddress("R4_bbc_mean",&(RR[0][4]));
  rdat->SetBranchAddress("R5_bbc_mean",&(RR[0][5]));
  rdat->SetBranchAddress("R6_bbc_mean",&(RR[0][6]));
  rdat->SetBranchAddress("R7_bbc_mean",&(RR[0][7]));
  rdat->SetBranchAddress("R8_bbc_mean",&(RR[0][8]));
  rdat->SetBranchAddress("R9_bbc_mean",&(RR[0][9]));
  rdat->SetBranchAddress("R1_zdc_mean",&(RR[1][1]));
  rdat->SetBranchAddress("R2_zdc_mean",&(RR[1][2]));
  rdat->SetBranchAddress("R3_zdc_mean",&(RR[1][3]));
  rdat->SetBranchAddress("R4_zdc_mean",&(RR[1][4]));
  rdat->SetBranchAddress("R5_zdc_mean",&(RR[1][5]));
  rdat->SetBranchAddress("R6_zdc_mean",&(RR[1][6]));
  rdat->SetBranchAddress("R7_zdc_mean",&(RR[1][7]));
  rdat->SetBranchAddress("R8_zdc_mean",&(RR[1][8]));
  rdat->SetBranchAddress("R9_zdc_mean",&(RR[1][9]));
  rdat->SetBranchAddress("R1_vpd_mean",&(RR[2][1]));
  rdat->SetBranchAddress("R2_vpd_mean",&(RR[2][2]));
  rdat->SetBranchAddress("R3_vpd_mean",&(RR[2][3]));
  rdat->SetBranchAddress("R4_vpd_mean",&(RR[2][4]));
  rdat->SetBranchAddress("R5_vpd_mean",&(RR[2][5]));
  rdat->SetBranchAddress("R6_vpd_mean",&(RR[2][6]));
  rdat->SetBranchAddress("R7_vpd_mean",&(RR[2][7]));
  rdat->SetBranchAddress("R8_vpd_mean",&(RR[2][8]));
  rdat->SetBranchAddress("R9_vpd_mean",&(RR[2][9]));

  rdat->SetBranchAddress("R1_bbc_mean_err",&(RR_err[0][1]));
  rdat->SetBranchAddress("R2_bbc_mean_err",&(RR_err[0][2]));
  rdat->SetBranchAddress("R3_bbc_mean_err",&(RR_err[0][3]));
  rdat->SetBranchAddress("R4_bbc_mean_err",&(RR_err[0][4]));
  rdat->SetBranchAddress("R5_bbc_mean_err",&(RR_err[0][5]));
  rdat->SetBranchAddress("R6_bbc_mean_err",&(RR_err[0][6]));
  rdat->SetBranchAddress("R7_bbc_mean_err",&(RR_err[0][7]));
  rdat->SetBranchAddress("R8_bbc_mean_err",&(RR_err[0][8]));
  rdat->SetBranchAddress("R9_bbc_mean_err",&(RR_err[0][9]));
  rdat->SetBranchAddress("R1_zdc_mean_err",&(RR_err[1][1]));
  rdat->SetBranchAddress("R2_zdc_mean_err",&(RR_err[1][2]));
  rdat->SetBranchAddress("R3_zdc_mean_err",&(RR_err[1][3]));
  rdat->SetBranchAddress("R4_zdc_mean_err",&(RR_err[1][4]));
  rdat->SetBranchAddress("R5_zdc_mean_err",&(RR_err[1][5]));
  rdat->SetBranchAddress("R6_zdc_mean_err",&(RR_err[1][6]));
  rdat->SetBranchAddress("R7_zdc_mean_err",&(RR_err[1][7]));
  rdat->SetBranchAddress("R8_zdc_mean_err",&(RR_err[1][8]));
  rdat->SetBranchAddress("R9_zdc_mean_err",&(RR_err[1][9]));
  rdat->SetBranchAddress("R1_vpd_mean_err",&(RR_err[2][1]));
  rdat->SetBranchAddress("R2_vpd_mean_err",&(RR_err[2][2]));
  rdat->SetBranchAddress("R3_vpd_mean_err",&(RR_err[2][3]));
  rdat->SetBranchAddress("R4_vpd_mean_err",&(RR_err[2][4]));
  rdat->SetBranchAddress("R5_vpd_mean_err",&(RR_err[2][5]));
  rdat->SetBranchAddress("R6_vpd_mean_err",&(RR_err[2][6]));
  rdat->SetBranchAddress("R7_vpd_mean_err",&(RR_err[2][7]));
  rdat->SetBranchAddress("R8_vpd_mean_err",&(RR_err[2][8]));
  rdat->SetBranchAddress("R9_vpd_mean_err",&(RR_err[2][9]));

  rdat->SetBranchAddress("d_vz_e",&(d_vz[0]));
  rdat->SetBranchAddress("d_vz_w",&(d_vz[1]));
  rdat->SetBranchAddress("d_vz_x",&(d_vz[2]));
  rdat->SetBranchAddress("d_ew_bbc",&(d_xx[0][0]));
  rdat->SetBranchAddress("d_ew_zdc",&(d_xx[0][1]));
  rdat->SetBranchAddress("d_ew_vpd",&(d_xx[0][2]));
  rdat->SetBranchAddress("d_ex_bbc",&(d_xx[1][0]));
  rdat->SetBranchAddress("d_ex_zdc",&(d_xx[1][1]));
  rdat->SetBranchAddress("d_ex_vpd",&(d_xx[1][2]));
  rdat->SetBranchAddress("d_wx_bbc",&(d_xx[2][0]));
  rdat->SetBranchAddress("d_wx_zdc",&(d_xx[2][1]));
  rdat->SetBranchAddress("d_wx_vpd",&(d_xx[2][2]));

  sums->SetBranchAddress("i",&i_sums);
  sums->SetBranchAddress("runnum",&runnum_sums);
  sums->SetBranchAddress("fi",&fi_sums);
  sums->SetBranchAddress("fill",&fill_sums);
  sums->SetBranchAddress("t",&t_sums);
  sums->SetBranchAddress("tau",&tau);
  sums->SetBranchAddress("num_runs",&num_runs);
  sums->SetBranchAddress("bbce",&bbce);
  sums->SetBranchAddress("bbcw",&bbcw);
  sums->SetBranchAddress("bbcx",&bbcx);
  sums->SetBranchAddress("zdce",&zdce);
  sums->SetBranchAddress("zdcw",&zdcw);
  sums->SetBranchAddress("zdcx",&zdcx);
  sums->SetBranchAddress("vpde",&vpde);
  sums->SetBranchAddress("vpdw",&vpdw);
  sums->SetBranchAddress("vpdx",&vpdx);
  sums->SetBranchAddress("tot_bx",&tot_bx);
  sums->SetBranchAddress("pattern",&pattern);

  
  TFile * outfile = new TFile("rtree.root","RECREATE");
  TTree * rellum = new TTree("rellum","rellum");
  rellum->Branch("i",&i_sums,"i/I"); // run index
  rellum->Branch("runnum",&runnum_sums,"runnum/I"); // run number
  rellum->Branch("fi",&fi_sums,"fi/I"); // fill index
  rellum->Branch("fill",&fill_sums,"fill/I"); // fill number
  rellum->Branch("t",&t_sums,"t/D"); // run time (seconds)
  rellum->Branch("tau",&tau,"tau/D"); // total bXings / bXing rate (nominally == run time)
  rellum->Branch("num_runs",&num_runs,"num_runs/I"); // total no. runs in a fill
  rellum->Branch("bbce",&bbce,"bbce/D"); // multiples & accidentals corrected scaler counts
  rellum->Branch("bbcw",&bbcw,"bbcw/D");
  rellum->Branch("bbcx",&bbcx,"bbcx/D");
  rellum->Branch("zdce",&zdce,"zdce/D");
  rellum->Branch("zdcw",&zdcw,"zdcw/D");
  rellum->Branch("zdcx",&zdcx,"zdcx/D");
  rellum->Branch("vpde",&vpde,"vpde/D");
  rellum->Branch("vpdw",&vpdw,"vpdw/D");
  rellum->Branch("vpdx",&vpdx,"vpdx/D");
  rellum->Branch("tot_bx",&tot_bx,"tot_bx/D"); // total no. bXings
  rellum->Branch("R1_bbce",&(RRe[0][1]),"R1_bbce/F"); // bbce relative luminosity
  rellum->Branch("R2_bbce",&(RRe[0][2]),"R2_bbce/F");
  rellum->Branch("R3_bbce",&(RRe[0][3]),"R3_bbce/F");
  rellum->Branch("R4_bbce",&(RRe[0][4]),"R4_bbce/F");
  rellum->Branch("R5_bbce",&(RRe[0][5]),"R5_bbce/F");
  rellum->Branch("R6_bbce",&(RRe[0][6]),"R6_bbce/F");
  rellum->Branch("R7_bbce",&(RRe[0][7]),"R7_bbce/F");
  rellum->Branch("R8_bbce",&(RRe[0][8]),"R8_bbce/F");
  rellum->Branch("R9_bbce",&(RRe[0][9]),"R9_bbce/F");
  rellum->Branch("R1_zdce",&(RRe[1][1]),"R1_zdce/F"); // zdce relative luminosity
  rellum->Branch("R2_zdce",&(RRe[1][2]),"R2_zdce/F");
  rellum->Branch("R3_zdce",&(RRe[1][3]),"R3_zdce/F");
  rellum->Branch("R4_zdce",&(RRe[1][4]),"R4_zdce/F");
  rellum->Branch("R5_zdce",&(RRe[1][5]),"R5_zdce/F");
  rellum->Branch("R6_zdce",&(RRe[1][6]),"R6_zdce/F");
  rellum->Branch("R7_zdce",&(RRe[1][7]),"R7_zdce/F");
  rellum->Branch("R8_zdce",&(RRe[1][8]),"R8_zdce/F");
  rellum->Branch("R9_zdce",&(RRe[1][9]),"R9_zdce/F");
  rellum->Branch("R1_vpde",&(RRe[2][1]),"R1_vpde/F"); // vpde relative luminosity
  rellum->Branch("R2_vpde",&(RRe[2][2]),"R2_vpde/F");
  rellum->Branch("R3_vpde",&(RRe[2][3]),"R3_vpde/F");
  rellum->Branch("R4_vpde",&(RRe[2][4]),"R4_vpde/F");
  rellum->Branch("R5_vpde",&(RRe[2][5]),"R5_vpde/F");
  rellum->Branch("R6_vpde",&(RRe[2][6]),"R6_vpde/F");
  rellum->Branch("R7_vpde",&(RRe[2][7]),"R7_vpde/F");
  rellum->Branch("R8_vpde",&(RRe[2][8]),"R8_vpde/F");
  rellum->Branch("R9_vpde",&(RRe[2][9]),"R9_vpde/F");
  rellum->Branch("R1_bbcw",&(RRw[0][1]),"R1_bbcw/F"); // bbcw relative luminosity
  rellum->Branch("R2_bbcw",&(RRw[0][2]),"R2_bbcw/F");
  rellum->Branch("R3_bbcw",&(RRw[0][3]),"R3_bbcw/F");
  rellum->Branch("R4_bbcw",&(RRw[0][4]),"R4_bbcw/F");
  rellum->Branch("R5_bbcw",&(RRw[0][5]),"R5_bbcw/F");
  rellum->Branch("R6_bbcw",&(RRw[0][6]),"R6_bbcw/F");
  rellum->Branch("R7_bbcw",&(RRw[0][7]),"R7_bbcw/F");
  rellum->Branch("R8_bbcw",&(RRw[0][8]),"R8_bbcw/F");
  rellum->Branch("R9_bbcw",&(RRw[0][9]),"R9_bbcw/F");
  rellum->Branch("R1_zdcw",&(RRw[1][1]),"R1_zdcw/F"); // zdcw relative luminosity
  rellum->Branch("R2_zdcw",&(RRw[1][2]),"R2_zdcw/F");
  rellum->Branch("R3_zdcw",&(RRw[1][3]),"R3_zdcw/F");
  rellum->Branch("R4_zdcw",&(RRw[1][4]),"R4_zdcw/F");
  rellum->Branch("R5_zdcw",&(RRw[1][5]),"R5_zdcw/F");
  rellum->Branch("R6_zdcw",&(RRw[1][6]),"R6_zdcw/F");
  rellum->Branch("R7_zdcw",&(RRw[1][7]),"R7_zdcw/F");
  rellum->Branch("R8_zdcw",&(RRw[1][8]),"R8_zdcw/F");
  rellum->Branch("R9_zdcw",&(RRw[1][9]),"R9_zdcw/F");
  rellum->Branch("R1_vpdw",&(RRw[2][1]),"R1_vpdw/F"); // vpdw relative luminosity
  rellum->Branch("R2_vpdw",&(RRw[2][2]),"R2_vpdw/F");
  rellum->Branch("R3_vpdw",&(RRw[2][3]),"R3_vpdw/F");
  rellum->Branch("R4_vpdw",&(RRw[2][4]),"R4_vpdw/F");
  rellum->Branch("R5_vpdw",&(RRw[2][5]),"R5_vpdw/F");
  rellum->Branch("R6_vpdw",&(RRw[2][6]),"R6_vpdw/F");
  rellum->Branch("R7_vpdw",&(RRw[2][7]),"R7_vpdw/F");
  rellum->Branch("R8_vpdw",&(RRw[2][8]),"R8_vpdw/F");
  rellum->Branch("R9_vpdw",&(RRw[2][9]),"R9_vpdw/F");
  rellum->Branch("R1_bbcx",&(RRx[0][1]),"R1_bbcx/F"); // bbcx relative luminosity
  rellum->Branch("R2_bbcx",&(RRx[0][2]),"R2_bbcx/F");
  rellum->Branch("R3_bbcx",&(RRx[0][3]),"R3_bbcx/F");
  rellum->Branch("R4_bbcx",&(RRx[0][4]),"R4_bbcx/F");
  rellum->Branch("R5_bbcx",&(RRx[0][5]),"R5_bbcx/F");
  rellum->Branch("R6_bbcx",&(RRx[0][6]),"R6_bbcx/F");
  rellum->Branch("R7_bbcx",&(RRx[0][7]),"R7_bbcx/F");
  rellum->Branch("R8_bbcx",&(RRx[0][8]),"R8_bbcx/F");
  rellum->Branch("R9_bbcx",&(RRx[0][9]),"R9_bbcx/F");
  rellum->Branch("R1_zdcx",&(RRx[1][1]),"R1_zdcx/F"); // zdcx relative luminosity
  rellum->Branch("R2_zdcx",&(RRx[1][2]),"R2_zdcx/F");
  rellum->Branch("R3_zdcx",&(RRx[1][3]),"R3_zdcx/F");
  rellum->Branch("R4_zdcx",&(RRx[1][4]),"R4_zdcx/F");
  rellum->Branch("R5_zdcx",&(RRx[1][5]),"R5_zdcx/F");
  rellum->Branch("R6_zdcx",&(RRx[1][6]),"R6_zdcx/F");
  rellum->Branch("R7_zdcx",&(RRx[1][7]),"R7_zdcx/F");
  rellum->Branch("R8_zdcx",&(RRx[1][8]),"R8_zdcx/F");
  rellum->Branch("R9_zdcx",&(RRx[1][9]),"R9_zdcx/F");
  rellum->Branch("R1_vpdx",&(RRx[2][1]),"R1_vpdx/F"); // vpdx relative luminosity
  rellum->Branch("R2_vpdx",&(RRx[2][2]),"R2_vpdx/F");
  rellum->Branch("R3_vpdx",&(RRx[2][3]),"R3_vpdx/F");
  rellum->Branch("R4_vpdx",&(RRx[2][4]),"R4_vpdx/F");
  rellum->Branch("R5_vpdx",&(RRx[2][5]),"R5_vpdx/F");
  rellum->Branch("R6_vpdx",&(RRx[2][6]),"R6_vpdx/F");
  rellum->Branch("R7_vpdx",&(RRx[2][7]),"R7_vpdx/F");
  rellum->Branch("R8_vpdx",&(RRx[2][8]),"R8_vpdx/F");
  rellum->Branch("R9_vpdx",&(RRx[2][9]),"R9_vpdx/F");
  rellum->Branch("R1_bbc_mean",&(RR[0][1]),"R1_bbc_mean/F"); // mean bbc relative luminosity
  rellum->Branch("R2_bbc_mean",&(RR[0][2]),"R2_bbc_mean/F");
  rellum->Branch("R3_bbc_mean",&(RR[0][3]),"R3_bbc_mean/F");
  rellum->Branch("R4_bbc_mean",&(RR[0][4]),"R4_bbc_mean/F");
  rellum->Branch("R5_bbc_mean",&(RR[0][5]),"R5_bbc_mean/F");
  rellum->Branch("R6_bbc_mean",&(RR[0][6]),"R6_bbc_mean/F");
  rellum->Branch("R7_bbc_mean",&(RR[0][7]),"R7_bbc_mean/F");
  rellum->Branch("R8_bbc_mean",&(RR[0][8]),"R8_bbc_mean/F");
  rellum->Branch("R9_bbc_mean",&(RR[0][9]),"R9_bbc_mean/F");
  rellum->Branch("R1_zdc_mean",&(RR[1][1]),"R1_zdc_mean/F"); // mean zdc relative luminosity
  rellum->Branch("R2_zdc_mean",&(RR[1][2]),"R2_zdc_mean/F");
  rellum->Branch("R3_zdc_mean",&(RR[1][3]),"R3_zdc_mean/F");
  rellum->Branch("R4_zdc_mean",&(RR[1][4]),"R4_zdc_mean/F");
  rellum->Branch("R5_zdc_mean",&(RR[1][5]),"R5_zdc_mean/F");
  rellum->Branch("R6_zdc_mean",&(RR[1][6]),"R6_zdc_mean/F");
  rellum->Branch("R7_zdc_mean",&(RR[1][7]),"R7_zdc_mean/F");
  rellum->Branch("R8_zdc_mean",&(RR[1][8]),"R8_zdc_mean/F");
  rellum->Branch("R9_zdc_mean",&(RR[1][9]),"R9_zdc_mean/F");
  rellum->Branch("R1_vpd_mean",&(RR[2][1]),"R1_vpd_mean/F"); // mean vpd relative luminosity
  rellum->Branch("R2_vpd_mean",&(RR[2][2]),"R2_vpd_mean/F");
  rellum->Branch("R3_vpd_mean",&(RR[2][3]),"R3_vpd_mean/F");
  rellum->Branch("R4_vpd_mean",&(RR[2][4]),"R4_vpd_mean/F");
  rellum->Branch("R5_vpd_mean",&(RR[2][5]),"R5_vpd_mean/F");
  rellum->Branch("R6_vpd_mean",&(RR[2][6]),"R6_vpd_mean/F");
  rellum->Branch("R7_vpd_mean",&(RR[2][7]),"R7_vpd_mean/F");
  rellum->Branch("R8_vpd_mean",&(RR[2][8]),"R8_vpd_mean/F");
  rellum->Branch("R9_vpd_mean",&(RR[2][9]),"R9_vpd_mean/F");
  rellum->Branch("R1_bbc_mean_err",&(RR_err[0][1]),"R1_bbc_mean_err/F"); // mean_err bbc rellum err
  rellum->Branch("R2_bbc_mean_err",&(RR_err[0][2]),"R2_bbc_mean_err/F");
  rellum->Branch("R3_bbc_mean_err",&(RR_err[0][3]),"R3_bbc_mean_err/F");
  rellum->Branch("R4_bbc_mean_err",&(RR_err[0][4]),"R4_bbc_mean_err/F");
  rellum->Branch("R5_bbc_mean_err",&(RR_err[0][5]),"R5_bbc_mean_err/F");
  rellum->Branch("R6_bbc_mean_err",&(RR_err[0][6]),"R6_bbc_mean_err/F");
  rellum->Branch("R7_bbc_mean_err",&(RR_err[0][7]),"R7_bbc_mean_err/F");
  rellum->Branch("R8_bbc_mean_err",&(RR_err[0][8]),"R8_bbc_mean_err/F");
  rellum->Branch("R9_bbc_mean_err",&(RR_err[0][9]),"R9_bbc_mean_err/F");
  rellum->Branch("R1_zdc_mean_err",&(RR_err[1][1]),"R1_zdc_mean_err/F"); // mean_err zdc rellum err
  rellum->Branch("R2_zdc_mean_err",&(RR_err[1][2]),"R2_zdc_mean_err/F");
  rellum->Branch("R3_zdc_mean_err",&(RR_err[1][3]),"R3_zdc_mean_err/F");
  rellum->Branch("R4_zdc_mean_err",&(RR_err[1][4]),"R4_zdc_mean_err/F");
  rellum->Branch("R5_zdc_mean_err",&(RR_err[1][5]),"R5_zdc_mean_err/F");
  rellum->Branch("R6_zdc_mean_err",&(RR_err[1][6]),"R6_zdc_mean_err/F");
  rellum->Branch("R7_zdc_mean_err",&(RR_err[1][7]),"R7_zdc_mean_err/F");
  rellum->Branch("R8_zdc_mean_err",&(RR_err[1][8]),"R8_zdc_mean_err/F");
  rellum->Branch("R9_zdc_mean_err",&(RR_err[1][9]),"R9_zdc_mean_err/F");
  rellum->Branch("R1_vpd_mean_err",&(RR_err[2][1]),"R1_vpd_mean_err/F"); // mean_err vpd rellum err
  rellum->Branch("R2_vpd_mean_err",&(RR_err[2][2]),"R2_vpd_mean_err/F");
  rellum->Branch("R3_vpd_mean_err",&(RR_err[2][3]),"R3_vpd_mean_err/F");
  rellum->Branch("R4_vpd_mean_err",&(RR_err[2][4]),"R4_vpd_mean_err/F");
  rellum->Branch("R5_vpd_mean_err",&(RR_err[2][5]),"R5_vpd_mean_err/F");
  rellum->Branch("R6_vpd_mean_err",&(RR_err[2][6]),"R6_vpd_mean_err/F");
  rellum->Branch("R7_vpd_mean_err",&(RR_err[2][7]),"R7_vpd_mean_err/F");
  rellum->Branch("R8_vpd_mean_err",&(RR_err[2][8]),"R8_vpd_mean_err/F");
  rellum->Branch("R9_vpd_mean_err",&(RR_err[2][9]),"R9_vpd_mean_err/F");
  rellum->Branch("d_vz_e",&(d_vz[0]),"d_vz_e/F"); // zdc-vpd diagnostic
  rellum->Branch("d_vz_w",&(d_vz[1]),"d_vz_w/F");
  rellum->Branch("d_vz_x",&(d_vz[2]),"d_vz_x/F");
  rellum->Branch("d_ew_bbc",&(d_xx[0][0]),"d_ew_bbc/F"); // east-west diagnostic
  rellum->Branch("d_ew_zdc",&(d_xx[0][1]),"d_ew_zdc/F");
  rellum->Branch("d_ew_vpd",&(d_xx[0][2]),"d_ew_vpd/F");
  rellum->Branch("d_ex_bbc",&(d_xx[1][0]),"d_ex_bbc/F"); // east-coin diagnostic
  rellum->Branch("d_ex_zdc",&(d_xx[1][1]),"d_ex_zdc/F");
  rellum->Branch("d_ex_vpd",&(d_xx[1][2]),"d_ex_vpd/F");
  rellum->Branch("d_wx_bbc",&(d_xx[2][0]),"d_wx_bbc/F"); // west-coin diagnostic
  rellum->Branch("d_wx_zdc",&(d_xx[2][1]),"d_wx_zdc/F");
  rellum->Branch("d_wx_vpd",&(d_xx[2][2]),"d_wx_vpd/F");
  rellum->Branch("isConsistent",&isConsistent,"isConsistent/O"); // true if diagnostics passed
  rellum->Branch("pattern",&pattern,"pattern/I"); // spin pattern no. (see sumTree.C)


  // diagnostic consistency bounds (see if statement in tree for loop below)
  // -- these bounds were determined by eye
  Float_t d_vz_cc[3];
  Float_t d_xx_cc[3][3];
  Float_t t_tau;
  d_vz_cc[0] = 0.006; // VPDE - ZDCE
  d_vz_cc[1] = 0.004; // VPDW - ZDCW
  d_vz_cc[2] = 0.006; // VPDX - ZDCX
  d_xx_cc[0][1] = 0.0015; // ZDCE - ZDCW
  d_xx_cc[0][2] = 0.002;  // VPDE - VPDW
  d_xx_cc[1][1] = 0.005;  // ZDCE - ZDCX
  d_xx_cc[1][2] = 0.002;  // VPDE - VPDX
  d_xx_cc[2][1] = 0.005;  // ZDCW - ZDCX
  d_xx_cc[2][2] = 0.002;  // VPDW - VPDX
  t_tau = 1.1;

  printf("\ndiagnostic cuts\n");
  printf("| VPDE - ZDCE | < %f\n",d_vz_cc[0]);
  printf("| VPDW - ZDCW | < %f\n",d_vz_cc[1]);
  printf("| VPDX - ZDCX | < %f\n",d_vz_cc[2]);
  printf("| ZDCE - ZDCW | < %f\n",d_xx_cc[0][1]);
  printf("| VPDE - VPDW | < %f\n",d_xx_cc[0][2]);
  printf("| ZDCE - ZDCX | < %f\n",d_xx_cc[1][1]);
  printf("| VPDE - VPDX | < %f\n",d_xx_cc[1][2]);
  printf("| ZDCW - ZDCX | < %f\n",d_xx_cc[2][1]);
  printf("| VPDW - VPDX | < %f\n",d_xx_cc[2][2]);
  printf("t / tau < %f\n",t_tau);

  // fill tree
  Int_t ent = rdat->GetEntries();
  for(Int_t e=0; e<ent; e++)
  {
    rdat->GetEntry(e);
    sums->GetEntry(e);

    // diagnostic pass check
    if( fabs(d_vz[0]) < d_vz_cc[0] &&
        fabs(d_vz[1]) < d_vz_cc[1] &&
        fabs(d_vz[2]) < d_vz_cc[2] &&
        fabs(d_xx[0][1]) < d_xx_cc[0][1] &&
        fabs(d_xx[0][2]) < d_xx_cc[0][2] &&
        fabs(d_xx[1][1]) < d_xx_cc[1][1] &&
        fabs(d_xx[1][2]) < d_xx_cc[1][2] &&
        fabs(d_xx[2][1]) < d_xx_cc[2][1] &&
        fabs(d_xx[2][2]) < d_xx_cc[2][2] &&
        t_sums/tau < t_tau )
    {
      isConsistent=1;
    }
    else isConsistent=0;


    // rdat and sum tree consistency check
    if( (i_rdat != i_sums) ||
        (runnum_rdat != runnum_sums) ||
        (fill_rdat != fill_sums) ||
        (t_rdat != t_sums))
    {
      fprintf(stderr,"ERROR: trees not consistent (either i,runnum,fill,t");
      return;
    }
    else rellum->Fill();
  };

  rellum->Write();

  printf("\nrellum tree written to rtree.root\n");
};
