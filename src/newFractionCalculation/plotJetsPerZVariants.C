// Plot the jets-per-Z ratio across the calculateRAA.C analysis switches.
//
// Expects one output file per configuration, produced by running calculateRAA.C
// with the switches set and copying the output aside. e.g.
//
//   root -l -b -q 'runVariant.C'   with
//     .L calculateRAA.C
//     doUnfolding=true; doFakeJetSubtraction=true; useMinBiasOnly=false;
//     calculateRAA();
//   cp rootFiles/JetsPerZ/histograms_JetsPerZ_lightJets_rebinned.root <inDir>/var_unf_fake.root
//
// Variants used here:
//   var_unf_fake     unfold + JEC fakes            (nominal)
//   var_unf_nofake   unfold, no subtraction
//   var_unf_nojec    unfold + non-JEC fakes
//   var_unf_pin      unfold + fakes pinned at 20 GeV
//   var_raw_fake     no unfolding + JEC fakes
//   var_raw_nofake   no unfolding, no subtraction
//   var_raw_pin      no unfolding + pinned fakes
//   var_mb_raw_fake  MinBias only, no unfolding, JEC fakes
//
// NOTE: all of these must come from the same N_iter. Iterative Bayes is only
// linear at N_iter = 1; at N >= 2 the prior is rebuilt from the data and the
// fake subtraction can raise the ratio, which is unphysical as a "correction".

const char *inDir  = "/tmp/claude-1000/-home-clayton-Analysis-code-bJetRaaAnalysis/0f4b0dfa-6f09-407c-a9e4-5f2eb201fe24/scratchpad/";
const char *outDir = "/home/clayton/Documents/nuclear/GroupMeeting/figures/2026-08-11/";

Int_t cc[5], cPin, cJec, cRaw;
const char *lbl[5]={"","PbPb 0-10%","PbPb 10-30%","PbPb 30-50%","PbPb 50-80%"};
int msty[5]={0,47,34,21,33};

void st(TH1 *h, Int_t col, int ms, int ls=1, double sz=0.9){
  h->SetLineColor(col); h->SetMarkerColor(col); h->SetLineWidth(2);
  h->SetLineStyle(ls); h->SetMarkerStyle(ms); h->SetMarkerSize(sz);
  h->SetStats(0); h->SetTitle("");
  h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.040); h->GetYaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitleOffset(1.40);
  h->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
}
void pd(TVirtualPad *p){ p->SetLeftMargin(0.16); p->SetBottomMargin(0.13);
                         p->SetTopMargin(0.06);  p->SetRightMargin(0.05); }
TH1D* G(const char *v, int ci, const char *nm){
  TFile *f=TFile::Open(Form("%svar_%s.root",inDir,v));
  if(!f||f->IsZombie()){ printf("cannot open var_%s.root\n",v); return nullptr; }
  TH1D *h=(TH1D*)f->Get(Form("r_C%d_fine",ci));
  if(!h){ printf("missing r_C%d_fine in var_%s\n",ci,v); return nullptr; }
  TH1D *o=(TH1D*)h->Clone(nm); o->SetDirectory(nullptr); f->Close(); return o;
}

// all four classes from one variant
void panel(const char *v, const char *tag, const char *out, double lo, double hi){
  TCanvas *c=new TCanvas(Form("c_%s",v),"",700,600); pd(c);
  TLegend *l=new TLegend(0.58,0.18,0.93,0.40);
  l->SetBorderSize(0); l->SetTextSize(0.036); l->SetFillStyle(0);
  for(int ci=1;ci<=4;ci++){
    TH1D *h=G(v,ci,Form("%s_p%d",v,ci)); if(!h) return;
    st(h,cc[ci],msty[ci]);
    h->GetYaxis()->SetTitle("jets per Z,  PbPb / pp");
    h->SetMinimum(lo); h->SetMaximum(hi);
    h->GetXaxis()->SetRangeUser(60,300);
    h->Draw(ci==1?"ep":"ep same");
    l->AddEntry(h,lbl[ci],"lp");
  }
  l->Draw();
  TLine *ln=new TLine(60,1,300,1); ln->SetLineStyle(2); ln->SetLineColor(kGray+1); ln->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.040); t.DrawLatex(0.20,0.88,tag);
  c->SaveAs(Form("%s%s",outDir,out));
}

void plotJetsPerZVariants(){
  gStyle->SetOptStat(0);
  cc[1]=TColor::GetColor("#D55E00"); cc[2]=TColor::GetColor("#0072B2");
  cc[3]=TColor::GetColor("#009E73"); cc[4]=TColor::GetColor("#E69F00");
  cJec=TColor::GetColor("#009E73"); cPin=TColor::GetColor("#CC79A7");
  cRaw=TColor::GetColor("#0072B2");

  panel("unf_nofake","No background subtraction","jetsPerZ_noBkgSub.pdf",   0.,1.25);
  panel("unf_fake",  "Fake jets subtracted",     "jetsPerZ_withBkgSub.pdf", 0.,1.25);
  panel("mb_raw_fake","MinBias only, no unfolding","jetsPerZ_minBiasOnly.pdf",0.,3.30);

  // ---- correction, all classes ----
  {
    TCanvas *c=new TCanvas("cEff","",700,600); pd(c);
    TLegend *l=new TLegend(0.58,0.18,0.93,0.40);
    l->SetBorderSize(0); l->SetTextSize(0.036); l->SetFillStyle(0);
    for(int ci=1;ci<=4;ci++){
      TH1D *a=G("unf_fake",ci,Form("ea%d",ci)), *b=G("unf_nofake",ci,Form("eb%d",ci));
      if(!a||!b) return;
      a->Divide(b); st(a,cc[ci],msty[ci]);
      a->GetYaxis()->SetTitle("(with bkg sub.) / (without)");
      a->SetMinimum(0.55); a->SetMaximum(1.08);
      a->GetXaxis()->SetRangeUser(60,300);
      a->Draw(ci==1?"ep":"ep same");
      l->AddEntry(a,lbl[ci],"lp");
    }
    l->Draw();
    TLine *ln=new TLine(60,1,300,1); ln->SetLineStyle(2); ln->SetLineColor(kGray+1); ln->Draw();
    c->SaveAs(Form("%sjetsPerZ_bkgSubEffect.pdf",outDir));
  }

  // ---- C4 alone, fakes on/off ----
  {
    TH1D *a=G("unf_fake",4,"c4a"), *b=G("unf_nofake",4,"c4b"); if(!a||!b) return;
    st(b,kGray+2,24); st(a,cc[4],20);
    TCanvas *c=new TCanvas("cC4","",700,600); pd(c);
    b->GetYaxis()->SetTitle("jets per Z, PbPb / pp   (50-80%)");
    b->SetMinimum(0.50); b->SetMaximum(0.90);
    b->GetXaxis()->SetRangeUser(60,300);
    b->Draw("ep"); a->Draw("ep same");
    TLegend *l=new TLegend(0.22,0.74,0.72,0.89);
    l->SetBorderSize(0); l->SetTextSize(0.038); l->SetFillStyle(0);
    l->AddEntry(b,"no background subtraction","lp");
    l->AddEntry(a,"fake jets subtracted","lp"); l->Draw();
    c->SaveAs(Form("%sjetsPerZ_C4.pdf",outDir));
  }

  // ---- factorisation: reco vs unfolded ----
  {
    TCanvas *c=new TCanvas("cFac","",700,600); pd(c);
    TLegend *l=new TLegend(0.50,0.18,0.94,0.42);
    l->SetBorderSize(0); l->SetTextSize(0.032); l->SetFillStyle(0);
    bool first=true;
    for(int ci : {1,4}){
      TH1D *r=G("raw_fake",ci,Form("fr%d",ci)), *rn=G("raw_nofake",ci,Form("frn%d",ci));
      TH1D *u=G("unf_fake",ci,Form("fu%d",ci)), *un=G("unf_nofake",ci,Form("fun%d",ci));
      if(!r||!rn||!u||!un) return;
      r->Divide(rn); u->Divide(un);
      st(r,cc[ci],msty[ci],2); st(u,cc[ci],msty[ci],1);
      r->GetYaxis()->SetTitle("(with bkg sub.) / (without)");
      r->SetMinimum(0.55); r->SetMaximum(1.10);
      r->GetXaxis()->SetRangeUser(60,300);
      r->Draw(first?"ep":"ep same"); first=false; u->Draw("ep same");
      l->AddEntry(r,Form("%s, reco level",lbl[ci]),"lp");
      l->AddEntry(u,Form("%s, unfolded",  lbl[ci]),"lp");
    }
    l->Draw();
    TLine *ln=new TLine(60,1,300,1); ln->SetLineStyle(2); ln->SetLineColor(kGray+1); ln->Draw();
    c->SaveAs(Form("%sjetsPerZ_factorise.pdf",outDir));
  }

  // ---- JEC / pin comparison, C1 and C4 ----
  {
    TCanvas *c=new TCanvas("cJec","",1200,560); c->Divide(2,1);
    int idx=0;
    for(int ci : {1,4}){
      idx++; c->cd(idx); pd(gPad);
      TH1D *o=G("unf_nofake",ci,Form("jo%d",ci)), *n=G("unf_nojec",ci,Form("jn%d",ci));
      TH1D *j=G("unf_fake",  ci,Form("jj%d",ci)), *p=G("unf_pin",  ci,Form("jp%d",ci));
      if(!o||!n||!j||!p) return;
      st(o,kGray+2,24); st(n,cRaw,25); st(j,cJec,20); st(p,cPin,22);
      o->GetYaxis()->SetTitle("jets per Z,  PbPb / pp");
      o->SetMinimum((ci==1)?0.10:0.45); o->SetMaximum((ci==1)?0.75:0.90);
      o->GetXaxis()->SetRangeUser(60,240);
      o->Draw("ep"); n->Draw("ep same"); j->Draw("ep same"); p->Draw("ep same");
      TLegend *l=new TLegend(0.44,0.62,0.95,0.88);
      l->SetBorderSize(0); l->SetTextSize(0.034); l->SetFillStyle(0);
      l->AddEntry(o,"no fake subtraction","lp");
      l->AddEntry(n,"fakes, no JEC","lp");
      l->AddEntry(j,"fakes, JEC","lp");
      l->AddEntry(p,"fakes pinned at 20 GeV","lp");
      l->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.048);
      t.DrawLatex(0.20,0.88, ci==1?"PbPb 0-10%":"PbPb 50-80%");
    }
    c->SaveAs(Form("%sjetsPerZ_JECeffect.pdf",outDir));
  }

  // ---- reco vs unfolded, incl. pinned ----
  {
    TCanvas *c=new TCanvas("cSign","",700,600); pd(c);
    TLegend *l=new TLegend(0.44,0.17,0.95,0.44);
    l->SetBorderSize(0); l->SetTextSize(0.030); l->SetFillStyle(0);
    bool first=true;
    for(int ci : {1,4}){
      TH1D *r=G("raw_fake",ci,Form("sr%d",ci)), *rn=G("raw_nofake",ci,Form("srn%d",ci));
      TH1D *u=G("unf_fake",ci,Form("su%d",ci)), *un=G("unf_nofake",ci,Form("sun%d",ci));
      if(!r||!rn||!u||!un) return;
      r->Divide(rn); u->Divide(un);
      st(r,cc[ci],msty[ci],2); st(u,cc[ci],msty[ci],1);
      r->GetYaxis()->SetTitle("(with fake sub.) / (without)");
      r->SetMinimum(0.30); r->SetMaximum(1.12);
      r->GetXaxis()->SetRangeUser(60,200);
      r->Draw(first?"ep":"ep same"); first=false; u->Draw("ep same");
      l->AddEntry(r,Form("%s, reco",lbl[ci]),"lp");
      l->AddEntry(u,Form("%s, unfolded",lbl[ci]),"lp");
      if(ci==1){
        TH1D *pr=G("raw_pin",ci,Form("spr%d",ci)), *prn=G("raw_nofake",ci,Form("sprn%d",ci));
        TH1D *pu=G("unf_pin",ci,Form("spu%d",ci)), *pun=G("unf_nofake",ci,Form("spun%d",ci));
        if(pr&&prn&&pu&&pun){
          pr->Divide(prn); pu->Divide(pun);
          st(pr,cPin,26,2); st(pu,cPin,22,1);
          pr->Draw("ep same"); pu->Draw("ep same");
          l->AddEntry(pr,"0-10%, pinned, reco","lp");
          l->AddEntry(pu,"0-10%, pinned, unfolded","lp");
        }
      }
    }
    l->Draw();
    TLine *ln=new TLine(60,1,200,1); ln->SetLineStyle(2); ln->SetLineColor(kGray+1); ln->Draw();
    c->SaveAs(Form("%sfakeCorrection_recoVsUnfolded.pdf",outDir));
  }
  printf("\ndone\n");
}
