// draws distribution of asymmetry determined from bunch fit
//
// requires the following files, all of which are produced from RunPatterns:
//  - fit_result.[numer].[denom].root
//  - pats/fit_result.[numer].[denom].pat*.root
//  - these files are then to be stored in the subdirectory 
//    output/omit_[Nomit]_bXings_after_aborts/

void DrawProjections(Int_t Nomit=1, const char * numer="zdce",
                                    const char * denom="vpdx")
{
  char filename[256];
  sprintf(filename,"output/omit_%d_bXings_after_aborts/fit_result.%s.%s.root",
    Nomit,numer,denom);
  TFile * infile = new TFile(filename,"READ");
  TH1D * hist_all = (TH1D*) infile->Get("/epsilon/epsi_a3_v_run");

  TFile * patfile[8];
  char patfile_n[8][256];
  Int_t pat_arr[8] = {13,14,23,24,31,32,41,42};
  Color_t colours[8] = {kOrange,kRed,kMagenta,kBlue,kCyan+1,kGreen+1,kYellow+2,kViolet-6}; 
  TH1D * hist_pat[8];
  for(Int_t pp=0; pp<8; pp++)
  {
    sprintf(patfile_n[pp],"output/omit_%d_bXings_after_aborts/pats/fit_result.%s.%s.pat%d.root",
      Nomit,numer,denom,pat_arr[pp]);
    patfile[pp] = new TFile(patfile_n[pp],"READ");
    hist_pat[pp] = (TH1D*) patfile[pp]->Get("/epsilon/epsi_a3_v_run");
  };


  char epsi_dist_all_n[256];
  char epsi_dist_all_t[256];
  sprintf(epsi_dist_all_n,"epsi_dist_all");
  sprintf(epsi_dist_all_t,"#varepsilon_{3} dist for %s/%s :: omitting %d bx after aborts",
    numer,denom,Nomit);
  TH1D * epsi_dist_all = new TH1D(epsi_dist_all_n,epsi_dist_all_t,50,-1e-3,1e-3);

  char epsi_dist_pat_n[8][256];
  char epsi_dist_pat_t[8][256];
  TH1D * epsi_dist_pat[8];
  for(Int_t pp=0; pp<8; pp++)
  {
    sprintf(epsi_dist_pat_n[pp],"epsi_dist_pat%d",pp);
    sprintf(epsi_dist_pat_t[pp],"#varepsilon_{3} dist for %s/%s :: omitting %d bx after aborts",
      numer,denom,Nomit);
    epsi_dist_pat[pp] = new TH1D(epsi_dist_pat_n[pp],epsi_dist_pat_t[pp],50,-1e-3,1e-3);
    epsi_dist_pat[pp]->SetLineColor(colours[pp]);
  };

  Double_t bc;

  for(Int_t b=1; b<=hist_all->GetNbinsX(); b++)
  {
    bc = hist_all->GetBinContent(b);
    epsi_dist_all->Fill(bc);
  };

  
  for(Int_t pp=0; pp<8; pp++)
  {
    for(Int_t b=1; b<=hist_pat[pp]->GetNbinsX(); b++)
    {
      bc = hist_pat[pp]->GetBinContent(b);
      if(bc!=0) epsi_dist_pat[pp]->Fill(bc);
    };
  };


  TCanvas * cc = new TCanvas("cc","cc",500,250);
  cc->SetGrid(1,0);
  gStyle->SetOptStat(1100);
  gStyle->SetStatFontSize(0.1);
  Float_t size = 0.05;
  epsi_dist_all->GetXaxis()->SetLabelSize(size);
  epsi_dist_all->GetYaxis()->SetLabelSize(size);
  epsi_dist_all->Draw();
  for(Int_t pp=0; pp<8; pp++)
  {
    epsi_dist_pat[pp]->GetXaxis()->SetLabelSize(size);
    epsi_dist_pat[pp]->GetYaxis()->SetLabelSize(size);
    epsi_dist_pat[pp]->Draw("same");
  };
  epsi_dist_all->Draw("same");

  char printname[64];
  sprintf(printname,"%d.png",Nomit);
  cc->Print(printname,"png");
    
};
