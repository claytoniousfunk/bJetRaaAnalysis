// Plot random-cone eta/phi UE maps from PbPb pfCandAnalyzer output.
//
// For each centrality bin produces:
//   1. 2D COLZ map: mean random-cone pT vs (eta, phi)
//   2. Eta projection: mean pT vs eta, averaged over phi
//   3. Phi projection: mean pT vs phi, averaged over eta
//
// Usage: root -l -b -q 'plotRandConeEtaPhi.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *inFile =
    "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
    "PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPseudoJets_pfCand_pseudoJetCandPtMin-0.0_2026-7-8_ultraFineCentBins.root";

const char *outDir =
    "../../../../figures/jetPt/pseudoJets/randConeEtaPhi/";

const int NCent = 17;
const char *centLabel[NCent] = { "0-80%", "0-5%", "5-10%", "10-15%", "15-20%",
				 "20-25%", "25-30%", "30-35%", "35-40%",
				 "40-45%", "45-50%", "50-55%", "55-60%",
				 "60-65%", "65-70%", "70-75%", "75-80%"};

const char *ratioLabel[NCent-2] = {"0-5%/5-10%","5-10%/10-15%","10-15%/15-20%",
				   "15-20%/20-25%","20-25%/25-30%","25-30%/30-35%",
				   "30-35%/35-40%","35-40%/40-45%","40-45%/45-50%",
				   "45-50%/50-55%","50-55%/55-60%","55-60%/60-65%",
				   "60-65%/65-70%","65-70%/70-75%","70-75%/75-80%"
};

// Okabe-Ito for projection overlays
const int projCols[NCent] = {
    kBlack,
    (Int_t)TColor::GetColor("#E69F00"),
    (Int_t)TColor::GetColor("#56B4E9"),
    (Int_t)TColor::GetColor("#009E73"),
    (Int_t)TColor::GetColor("#0072B2")
};

static void styleAxes2D(TProfile2D *h, const char *ztitle)
{
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle("#phi");
    h->GetZaxis()->SetTitle(ztitle);
    h->GetXaxis()->SetTitleSize(0.052); h->GetXaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetTitleSize(0.052); h->GetYaxis()->SetLabelSize(0.042);
    h->GetZaxis()->SetTitleSize(0.046); h->GetZaxis()->SetLabelSize(0.038);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetZaxis()->SetTitleOffset(1.4);
}

static void styleAxes1D(TH1 *h, const char *xtitle, const char *ytitle)
{
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitleSize(0.052); h->GetXaxis()->SetLabelSize(0.042);
    h->GetYaxis()->SetTitleSize(0.052); h->GetYaxis()->SetLabelSize(0.042);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.4);
}

void plotRandConeEtaPhi()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(64);

    TFile *f = TFile::Open(inFile);
    if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inFile); return; }

    gSystem->Exec(Form("mkdir -p %s", outDir));

    TProfile2D *hMap[NCent];
    TH1D       *hEta[NCent];
    TH1D       *hPhi[NCent];

    for(int ci = 0; ci < NCent; ci++){
        f->GetObject(Form("h_randConeEtaPhi_C%d", ci), hMap[ci]);
        if(!hMap[ci]){ printf("ERROR: missing h_randConeEtaPhi_C%d\n", ci); return; }
        hMap[ci] = (TProfile2D*) hMap[ci]->Clone(Form("hMap_C%d", ci));
        hMap[ci]->SetDirectory(nullptr);

        // ProfileX/Y average the mean pT across the collapsed axis,
        // giving the correct mean RC pT vs eta (or phi).
        // ProjectionX/Y would instead SUM cell means → unphysically large values.
        hEta[ci] = (TH1D*) hMap[ci]->ProfileX(Form("hEta_C%d", ci));
        hPhi[ci] = (TH1D*) hMap[ci]->ProfileY(Form("hPhi_C%d", ci));
        hEta[ci]->SetDirectory(nullptr);
        hPhi[ci]->SetDirectory(nullptr);
    }
    f->Close();

    // ------------------------------------------------------------------
    // 1. 2D COLZ maps — one canvas per centrality
    // ------------------------------------------------------------------
    const char *ztitle = "Mean random-cone p_{T} [GeV]";
    for(int ci = 0; ci < NCent; ci++){
        TCanvas *c = new TCanvas(Form("cMap_C%d", ci), "", 750, 600);
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.13);
        c->SetTopMargin(0.10);

        styleAxes2D(hMap[ci], ztitle);
        hMap[ci]->Draw("COLZ");

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.042);
        lat.DrawLatex(0.13, 0.925,
            Form("PbPb 5.02 TeV  random cone R=0.4  cent. %s", centLabel[ci]));

        TString out = TString(outDir) + Form("randConeEtaPhi_C%d.pdf", ci);
        c->SaveAs(out);
        printf("Saved %s\n", out.Data());
        delete c;
    }

    // ------------------------------------------------------------------
    // 2. Eta projections overlaid
    // ------------------------------------------------------------------
    {
        // double ylo = 1e30, yhi = -1e30;
        // for(int ci = 0; ci < NCent; ci++){
        //     for(int b = 1; b <= hEta[ci]->GetNbinsX(); b++){
        //         double v = hEta[ci]->GetBinContent(b);
        //         if(v > 0 && v < ylo) ylo = v;
        //         if(v > yhi) yhi = v;
        //     }
        // }
        // if(yhi <= 0) yhi = 1.;
        // if(ylo > 1e29) ylo = 0.;

        TCanvas *c = new TCanvas("cEta", "", 750, 580);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.10);

        for(int ci = 0; ci < NCent; ci++){
            styleAxes1D(hEta[ci], "#eta", ztitle);
            //hEta[ci]->SetLineColor(projCols[ci]);
	    //ci == 0 ? hEta[ci]->SetLineColor(kBlack) : 
            hEta[ci]->SetLineWidth(2);
            // hEta[ci]->GetYaxis()->SetRangeUser(ylo * 0.85, yhi * 1.15);
	    hEta[ci]->GetYaxis()->SetRangeUser(0,170);
            hEta[ci]->Draw(ci == 0 ? "hist e PLC" : "hist e same PLC");
        }

        TLegend *lg = new TLegend(0.36, 0.62, 0.74, 0.88);
        lg->SetBorderSize(0); lg->SetTextSize(0.03);
	lg->SetNColumns(3);
        for(int ci = 0; ci < NCent; ci++){
            lg->AddEntry(hPhi[ci], centLabel[ci], "l");
	    lg->Draw();
	}

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
        lat.DrawLatex(0.16, 0.925, "PbPb 5.02 TeV  random cone R=0.4  #eta projection");

        TString out = TString(outDir) + "randConeEtaPhi_etaProjection.pdf";
        c->SaveAs(out);
        printf("Saved %s\n", out.Data());
        delete c;
    }

    // ------------------------------------------------------------------
    // 3. Phi projections overlaid
    // ------------------------------------------------------------------
    {
        double ylo = 1e30, yhi = -1e30;
        for(int ci = 0; ci < NCent; ci++){
            for(int b = 1; b <= hPhi[ci]->GetNbinsX(); b++){
                double v = hPhi[ci]->GetBinContent(b);
                if(v > 0 && v < ylo) ylo = v;
                if(v > yhi) yhi = v;
            }
        }
        if(yhi <= 0) yhi = 1.;
        if(ylo > 1e29) ylo = 0.;

        TCanvas *c = new TCanvas("cPhi", "", 750, 580);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.10);

        for(int ci = 0; ci < NCent; ci++){
            styleAxes1D(hPhi[ci], "#phi", ztitle);
            //hPhi[ci]->SetLineColor(projCols[ci]);
            hPhi[ci]->SetLineWidth(2);
            //hPhi[ci]->GetYaxis()->SetRangeUser(ylo * 0.85, yhi * 1.15);
	    hPhi[ci]->GetYaxis()->SetRangeUser(0,170);
            hPhi[ci]->Draw(ci == 0 ? "hist e PLC" : "hist e same PLC");
        }

	TLegend *lg = new TLegend(0.36, 0.62, 0.74, 0.88);
        lg->SetBorderSize(0); lg->SetTextSize(0.03);
	lg->SetNColumns(3);
        for(int ci = 0; ci < NCent; ci++){
            lg->AddEntry(hEta[ci], centLabel[ci], "l");
	    lg->Draw();
	}

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
        lat.DrawLatex(0.16, 0.925, "PbPb 5.02 TeV  random cone R=0.4  #phi projection");

        TString out = TString(outDir) + "randConeEtaPhi_phiProjection.pdf";
        c->SaveAs(out);
        printf("Saved %s\n", out.Data());
        delete c;
    }

    // Eta centrality evolution ratios
    // ratio of cent(i+1) / cent(i) 

    {
      std::vector<TH1D*> ratios;
      for(int ci = 1; ci < NCent-1; ci++){
	TH1D* r = (TH1D*) hMap[ci]->ProjectionX()->Clone(Form("r_C%i",ci));
	r->Divide(hMap[ci]->ProjectionX(),hMap[ci+1]->ProjectionX(),1,1,"");
	ratios.push_back(r);
      }

      TCanvas *c = new TCanvas("etaRatios", "", 750, 580);
      c->SetLeftMargin(0.15);
      c->SetRightMargin(0.05);
      c->SetBottomMargin(0.14);
      c->SetTopMargin(0.10);
      c->cd();

      for(int ci = 1; ci < NCent-1; ci++){
	TH1D *r = (TH1D*) ratios.at(ci-1);
	r->SetLineWidth(2);
	r->GetYaxis()->SetRangeUser(0,2);
	r->GetYaxis()->SetTitle("<p_{T}^{RC}>(cent i) / <p_{T}^{RC}>(cent i+1)");
	r->GetYaxis()->SetTitleSize(0.04);
	r->GetXaxis()->SetTitle("#eta");
	r->Draw(ci == 1 ? "hist e PLC" : "hist e same PLC");
	
      }

      gPad->Update();

      TLegend *lg = new TLegend(0.16, 0.2, 0.74, 0.45);
      lg->SetBorderSize(0); lg->SetTextSize(0.03);
      lg->SetNColumns(3);
      for(int ci = 0; ci < NCent-2; ci++){
	lg->AddEntry(ratios.at(ci), ratioLabel[ci], "l");
      }

      lg->Draw();

      TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
      lat.DrawLatex(0.16, 0.925, "#eta projection centrality evolution");

      TString out = TString(outDir) + "randConeEtaPhi_etaProjection_centEvolution.pdf";
      c->SaveAs(out);
      printf("Saved %s\n", out.Data());
      delete c;



    }


    // Phi centrality evolution ratios
    // ratio of cent(i+1) / cent(i) 

    {
      std::vector<TH1D*> ratios;
      for(int ci = 1; ci < NCent-1; ci++){
	TH1D* r = (TH1D*) hMap[ci]->ProjectionY()->Clone(Form("r_C%i",ci));
	r->Divide(hMap[ci]->ProjectionY(),hMap[ci+1]->ProjectionY(),1,1,"");
	ratios.push_back(r);
      }

      TCanvas *c = new TCanvas("etaRatios", "", 750, 580);
      c->SetLeftMargin(0.15);
      c->SetRightMargin(0.05);
      c->SetBottomMargin(0.14);
      c->SetTopMargin(0.10);
      c->cd();
      
      for(int ci = 1; ci < NCent-1; ci++){
	TH1D *r = (TH1D*) ratios.at(ci-1);
	r->SetLineWidth(2);
	r->GetYaxis()->SetRangeUser(0,2);
	r->GetYaxis()->SetTitle("<p_{T}^{RC}>(cent i) / <p_{T}^{RC}>(cent i+1)");
	r->GetYaxis()->SetTitleSize(0.04);
	r->GetXaxis()->SetTitle("#phi");
	r->Draw(ci == 1 ? "hist e PLC" : "hist e same PLC");
      }

      TLegend *lg = new TLegend(0.16, 0.2, 0.84, 0.45);
      lg->SetBorderSize(0); lg->SetTextSize(0.03);
      lg->SetNColumns(3);
      for(int ci = 0; ci < NCent-2; ci++){
	lg->AddEntry(ratios.at(ci), ratioLabel[ci], "l");
	lg->Draw();
      }

      TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
      lat.DrawLatex(0.16, 0.925, "#phi projection centrality evolution");

      TString out = TString(outDir) + "randConeEtaPhi_phiProjection_centEvolution.pdf";
      c->SaveAs(out);
      printf("Saved %s\n", out.Data());
      delete c;


    }

    

    for(int ci = 0; ci < NCent; ci++){
        delete hMap[ci]; delete hEta[ci]; delete hPhi[ci];
    }

    printf("\nDone. Output in %s\n", outDir);
}
