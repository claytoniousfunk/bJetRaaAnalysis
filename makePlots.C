void makePlots(){
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    TFile *f15  = new TFile("rootFiles/scanningOuput/pp/pp_HighEGJet_Jet15HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-16.root");
    TFile *f30  = new TFile("rootFiles/scanningOuput/pp/pp_HighEGJet_Jet30HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-16.root");
    TFile *f100 = new TFile("rootFiles/scanningOuput/pp/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-6-16.root");

    const char* outDir = "/home/clayton/Documents/nuclear/claudeOutput/jet_comparison_plots/";
    int c15=kBlue+1, c30=kGreen+2, c100=kRed+1;

    // histogram names, x-axis labels, output file, logy, normalise
    struct PlotDef { const char* hname; const char* xtitle; const char* fname; bool logy; bool norm; };
    PlotDef plots[] = {
        {"h_inclRecoJetPt",              "Jet p_{T} [GeV]",              "inclJetPt.pdf",     true,  true},
        {"h_muTaggedRecoJetPt_triggerOn","#mu-tagged jet p_{T} [GeV]",   "muTaggedJetPt.pdf", true,  true},
        {"h_inclRecoJetEta",             "Jet #eta",                     "inclJetEta.pdf",    false, true},
        {"h_inclRecoJetPhi",             "Jet #phi",                     "inclJetPhi.pdf",    false, true},
        {"h_inclMuPt",                   "Muon p_{T} [GeV]",             "inclMuPt.pdf",      true,  true},
        // muPt and muPtrel handled separately via ProjectionX() below

        // muJetDR handled separately via ProjectionX() below
        {"h_dimuonMass",                 "m_{#mu#mu} [GeV]",             "dimuonMass.pdf",    false, false},
    };

    for(auto& p : plots){
        TCanvas* c = new TCanvas("c","",800,600);
        if(p.logy) c->SetLogy();

        TH1* h15  = (TH1*)f15 ->Get(p.hname); h15 ->Sumw2();
        TH1* h30  = (TH1*)f30 ->Get(p.hname); h30 ->Sumw2();
        TH1* h100 = (TH1*)f100->Get(p.hname); h100->Sumw2();

        for(TH1* h : {h15, h30, h100}){
            h->SetTitle("");
            h->GetXaxis()->SetTitle(p.xtitle);
            h->GetYaxis()->SetTitle(p.norm ? "Normalised entries / bin" : "Entries / bin");
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetLabelSize(0.04); h->GetYaxis()->SetLabelSize(0.04);
        }
        h15 ->SetLineColor(c15);  h15 ->SetLineWidth(2);
        h30 ->SetLineColor(c30);  h30 ->SetLineWidth(2);
        h100->SetLineColor(c100); h100->SetLineWidth(2);

        if(p.norm){
            if(h15 ->Integral()>0) h15 ->Scale(1./h15 ->Integral());
            if(h30 ->Integral()>0) h30 ->Scale(1./h30 ->Integral());
            if(h100->Integral()>0) h100->Scale(1./h100->Integral());
        }

        double ymax = TMath::Max(TMath::Max(h15->GetMaximum(), h30->GetMaximum()), h100->GetMaximum());
        if(!p.logy) h15->SetMaximum(ymax*1.4);
        h15->Draw("hist"); h30->Draw("hist same"); h100->Draw("hist same");

        TLegend* lg = new TLegend(0.55,0.65,0.88,0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(h15, "Jet15","l");
        lg->AddEntry(h30, "Jet30","l");
        lg->AddEntry(h100,"Jet100","l");
        lg->Draw();

        c->SaveAs(Form("%s%s", outDir, p.fname));
        delete lg; delete c;
    }

    // --- muon pT and pTrel: ProjectionX of TH2 (X=mu var, Y=jet pT) ---
    auto drawProj = [&](const char* hname, const char* xtitle, const char* fname, bool logy){
        TCanvas* c = new TCanvas("c","",800,600);
        if(logy) c->SetLogy();
        TH1* h15  = ((TH2*)f15 ->Get(hname))->ProjectionX(Form("%s_px15" ,hname));
        TH1* h30  = ((TH2*)f30 ->Get(hname))->ProjectionX(Form("%s_px30" ,hname));
        TH1* h100 = ((TH2*)f100->Get(hname))->ProjectionX(Form("%s_px100",hname));
        for(TH1* h : {h15, h30, h100}){
            h->SetTitle("");
            h->GetXaxis()->SetTitle(xtitle);
            h->GetYaxis()->SetTitle("Normalised entries / bin");
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetLabelSize(0.04); h->GetYaxis()->SetLabelSize(0.04);
        }
        h15 ->SetLineColor(c15);  h15 ->SetLineWidth(2);
        h30 ->SetLineColor(c30);  h30 ->SetLineWidth(2);
        h100->SetLineColor(c100); h100->SetLineWidth(2);
        if(h15 ->Integral()>0) h15 ->Scale(1./h15 ->Integral());
        if(h30 ->Integral()>0) h30 ->Scale(1./h30 ->Integral());
        if(h100->Integral()>0) h100->Scale(1./h100->Integral());
        double ymax = TMath::Max(TMath::Max(h15->GetMaximum(),h30->GetMaximum()),h100->GetMaximum());
        if(!logy) h15->SetMaximum(ymax*1.4);
        h15->Draw("hist"); h30->Draw("hist same"); h100->Draw("hist same");
        TLegend* lg = new TLegend(0.55,0.65,0.88,0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(h15,"Jet15","l"); lg->AddEntry(h30,"Jet30","l"); lg->AddEntry(h100,"Jet100","l");
        lg->Draw();
        c->SaveAs(Form("%s%s",outDir,fname));
        delete lg; delete c;
    };

    drawProj("h_muJetDr_recoJetPt",
             "#DeltaR(#mu, jet)",                         "muJetDR.pdf",     false);
    drawProj("h_mupt_recoJetPt_inclRecoMuonTag_triggerOn",
             "Muon p_{T} [GeV] (#mu-tagged, trig on)", "muPt_trigOn.pdf", true);
    drawProj("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn",
             "Muon p_{T}^{rel} [GeV]",                 "muPtRel.pdf",     true);
}
