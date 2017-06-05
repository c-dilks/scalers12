// plots number of SAME spin minus number OPPOSITE spin for 'num_ap_bxings' after
// each abort gap, in order to investigate biases in the low after-pulsing region (which
// is thought to occur right after each abort gap)


void SOasym(Int_t num_ap_bxings=4) {
  // open root file
  TFile * infile = new TFile("counts.root","READ");
  TTree * tr = (TTree*) infile->Get("sca");

  // open tree
  Int_t i,bx,blue,yell;
  Bool_t kicked;
  tr->SetBranchAddress("i",&i);
  tr->SetBranchAddress("bx",&bx);
  tr->SetBranchAddress("blue",&blue);
  tr->SetBranchAddress("yell",&yell);
  tr->SetBranchAddress("kicked",&kicked);

  // initialise some vars
  const Int_t NRUNS = 2000; // assumed
  Int_t asym[NRUNS][2]; // [run idx] [0=after_first_gap 1=after_second_gap]
  const Int_t afpulsebx[2] = {0,40};

  for(int d=0; d<NRUNS; d++) {
    for(int a=0; a<2; a++) {
      asym[d][a] = -1000;
    };
  };

  // loop through tree, filling up 'asym' var as we go
  Int_t runcount=0;
  for(int x=0; x<tr->GetEntries(); x++) {
    tr->GetEntry(x);
    if(i>=0 && i<NRUNS) {
      if(!kicked) {

        for(int a=0; a<2; a++) {
          //printf("i=%d a=%d bx=%d\n",i,a,bx);
          if(asym[i][a]==-1000) asym[i][a]=0;

          if(bx >= afpulsebx[a] && 
             bx <  afpulsebx[a]+num_ap_bxings) {
            asym[i][a] += blue*yell;
            //printf("bx=%d blue=%d yell=%d i=%d a=%d asym=%d\n",
              //bx,blue,yell,i,a,asym[i][a]);
          };
        };
      };
    } else {
      fprintf(stderr,"run index i=%d out of range!!!\n",i);
      return;
    };
  };

  // define graphs
  TGraph * gr[2];
  TGraph * gr_both;
  TString gr_n[2];
  TString gr_t[2];
  TString gr_both_t;
  Int_t gr_cnt[2];
  Int_t gr_both_cnt;
  for(int a=0; a<2; a++) {
    gr[a] = new TGraph();
    gr_n[a] = Form("asym_vs_run_%s",a==0?"L":"H");
    gr_t[a] = Form("%s post-abort-gap %d-bXings S/O asym vs. run index",
                   a==0?"1st":"2nd",num_ap_bxings);
    gr[a]->SetName(gr_n[a].Data());
    gr[a]->SetTitle(gr_t[a].Data());
    gr[a]->SetMarkerStyle(kFullCircle);
    gr_cnt[a]=0;
  };
  gr_both = new TGraph();
  gr_both->SetName("asym_vs_run_for_both_gaps");
  gr_both_t = Form("both post-abort-gap %d-bXings S/O asym vs. run index",num_ap_bxings);
  gr_both->SetTitle(gr_both_t.Data());
  gr_both_cnt = 0;
  gr_both->SetMarkerStyle(kFullCircle);

  // fill graphs
  for(int d=0; d<NRUNS; d++) {
    for(int a=0; a<2; a++) {
      if(asym[d][a]>-1000) {
        //printf("asym[%d][%d] = %d\n",d,a,asym[d][a]);
        gr[a]->SetPoint(gr_cnt[a]++,d,asym[d][a]);
      };
    };
    if(asym[d][0]>-1000 && asym[d][1]>-1000)
      gr_both->SetPoint(gr_both_cnt++,d,asym[d][0]+asym[d][1]);
  };

  // draw graphs
  TCanvas * canv = new TCanvas("canv","canv",1000,1000);
  canv->Divide(1,3);
  for(int a=0; a<3; a++) {
    canv->cd(a+1);
    canv->GetPad(a+1)->SetGrid(1,1);
    if(a<2) {
      gr[a]->GetYaxis()->SetRangeUser(-1*(num_ap_bxings+1),num_ap_bxings+1);
      gr[a]->Draw("APE");
    } else {
      gr_both->GetYaxis()->SetRangeUser(-1*(2*num_ap_bxings+1),2*num_ap_bxings+1);
      gr_both->Draw("APE");
    };
  };
};
