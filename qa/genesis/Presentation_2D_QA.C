// Presentation-quality standalone versions of the two 2D QA plots from QA_JewelTrees.C
// (truth-vs-reco leading jet pT correlation, average calo eta-phi energy image),
// for two JEWEL productions side by side. Larger canvas, bigger fonts, PNG + PDF output.
//
// Usage:
//   root -l -b -q 'Presentation_2D_QA.C("<outdir>", "<spec>", "<energy label>")'
// with spec = "name|legend label|tree file glob;name|legend label|tree file glob".
// Defaults reproduce the v1 vacuum/norecoil plots.

#include "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.572/rootmacros/sPhenixStyle.C"

#include <TCanvas.h>
#include <TChain.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TStyle.h>

#include <algorithm>
#include <string>
#include <vector>

struct Sample
{
  std::string name;
  std::string label;
  TChain *towertree;
  TChain *jettree;
};

// Three header lines in the pad's (enlarged) top margin, above the frame, so they can never
// cover data -- a 2D map has no guaranteed-empty region inside the frame the way a falling
// 1D spectrum does:
//   sample label
//   sPHENIX Simulation
//   <energy label>
void DrawHeader(const std::string &sampleLabel, const std::string &energyLabel)
{
  const double x = gPad->GetLeftMargin();
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.048);
  l.DrawLatex(x, 0.935, sampleLabel.c_str());
  l.SetTextFont(62);
  l.SetTextSize(0.042);
  l.DrawLatex(x, 0.880, "#it{#bf{sPHENIX}} Simulation");
  l.SetTextFont(42);
  l.SetTextSize(0.038);
  l.DrawLatex(x, 0.828, energyLabel.c_str());
}

void Presentation_2D_QA(
    const std::string &outdir = "/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/qa_plots",
    const std::string &spec = "vacuum|JEWEL vacuum|/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/vacuum/trees/*.root;"
                              "norecoil|JEWEL no-recoil|/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/norecoil/trees/*.root",
    const std::string &kEnergyLabel = "JEWEL p+p, #sqrt{s} = 200 GeV")
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kBird);
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gSystem->mkdir(outdir.c_str(), kTRUE);

  std::vector<Sample> samples;
  TObjArray *entries = TString(spec).Tokenize(";");
  for (int i = 0; i < entries->GetEntries(); i++)
  {
    TObjArray *f = ((TObjString *) entries->At(i))->GetString().Tokenize("|");
    if (f->GetEntries() != 3) continue;
    Sample s;
    s.name = ((TObjString *) f->At(0))->GetString().Data();
    s.label = ((TObjString *) f->At(1))->GetString().Data();
    s.towertree = new TChain("towertree");
    s.towertree->Add(((TObjString *) f->At(2))->GetString().Data());
    s.jettree = new TChain("jettree");
    s.jettree->Add(((TObjString *) f->At(2))->GetString().Data());
    samples.push_back(s);
  }

  //======================================================================
  // Truth-vs-reco leading jet pT correlation, side-by-side, big canvas
  //======================================================================
  {
    TCanvas *c1 = new TCanvas("c_truthvsreco", "c_truthvsreco", 1700, 950);
    c1->Divide(2, 1);
    for (size_t i = 0; i < samples.size(); i++)
    {
      c1->cd(i + 1);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.20);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.24);
      auto &s = samples[i];
      TH2D *h2 = new TH2D(Form("h2_truthvsreco_%s", s.name.c_str()),
                           ";leading truth jet p_{T} [GeV];leading reco jet p_{T} [GeV]",
                           25, 0, 60, 25, 0, 60);
      s.jettree->Draw(Form("Max$(rjet_pt):Max$(tjet_pt)>>%s", h2->GetName()),
                       "Length$(tjet_pt)>0 && Length$(rjet_pt)>0", "goff colz");
      h2->GetXaxis()->SetTitleOffset(1.2);
      h2->GetYaxis()->SetTitleOffset(1.3);
      h2->GetZaxis()->SetTitleOffset(1.5);
      h2->GetZaxis()->SetLabelSize(0.04);
      h2->Draw("colz");
      TLine line;
      line.SetLineStyle(2);
      line.SetLineWidth(2);
      line.DrawLine(0, 0, 60, 60);
      DrawHeader(s.label, kEnergyLabel);
    }
    c1->SaveAs((outdir + "/presentation_truth_vs_reco_pt.png").c_str());
    c1->SaveAs((outdir + "/presentation_truth_vs_reco_pt.pdf").c_str());
  }

  //======================================================================
  // Average calo image (eta-phi energy map), side-by-side, big canvas
  //======================================================================
  {
    TCanvas *c2 = new TCanvas("c_avgimg", "c_avgimg", 1700, 950);
    c2->Divide(2, 1);
    double pi = TMath::Pi();
    for (size_t i = 0; i < samples.size(); i++)
    {
      c2->cd(i + 1);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.20);
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.24);
      auto &s = samples[i];
      TH2D *h2evt = new TH2D(Form("h2_evt_%s", s.name.c_str()), ";#eta;#phi [rad]",
                              24, -1.1, 1.1, 64, 0, 2 * pi);

      std::vector<float> *em = nullptr, *ih = nullptr, *oh = nullptr;
      std::vector<int> *etabin = nullptr, *phibin = nullptr;
      s.towertree->SetBranchAddress("emcalenergy", &em);
      s.towertree->SetBranchAddress("ihcalenergy", &ih);
      s.towertree->SetBranchAddress("ohcalenergy", &oh);
      s.towertree->SetBranchAddress("etabin", &etabin);
      s.towertree->SetBranchAddress("phibin", &phibin);

      // Pick a single, representative event: the one with the highest total tower energy
      // (so the display isn't a near-empty low-activity event).
      Long64_t n = s.towertree->GetEntries();
      Long64_t bestEntry = 0;
      double bestTot = -1;
      for (Long64_t e = 0; e < n; e++)
      {
        s.towertree->GetEntry(e);
        double tot = 0;
        for (size_t k = 0; k < em->size(); k++)
          tot += std::max(0.f, (*em)[k]) + std::max(0.f, (*ih)[k]) + std::max(0.f, (*oh)[k]);
        if (tot > bestTot)
        {
          bestTot = tot;
          bestEntry = e;
        }
      }

      s.towertree->GetEntry(bestEntry);
      for (size_t k = 0; k < em->size(); k++)
      {
        double tot = std::max(0.f, (*em)[k]) + std::max(0.f, (*ih)[k]) + std::max(0.f, (*oh)[k]);
        h2evt->SetBinContent((*etabin)[k] + 1, (*phibin)[k] + 1,
                              h2evt->GetBinContent((*etabin)[k] + 1, (*phibin)[k] + 1) + tot);
      }
      h2evt->GetZaxis()->SetTitle("E [GeV]");
      h2evt->GetXaxis()->SetTitleOffset(1.2);
      h2evt->GetYaxis()->SetTitleOffset(1.2);
      h2evt->GetZaxis()->SetTitleOffset(1.5);
      h2evt->GetZaxis()->SetLabelSize(0.04);
      h2evt->Draw("colz");
      DrawHeader(s.label, kEnergyLabel);
      s.towertree->ResetBranchAddresses();
    }
    c2->SaveAs((outdir + "/presentation_single_event_calo_image.png").c_str());
    c2->SaveAs((outdir + "/presentation_single_event_calo_image.pdf").c_str());
  }

  std::cout << "Wrote presentation_truth_vs_reco_pt.{png,pdf} and presentation_single_event_calo_image.{png,pdf} to "
            << outdir << std::endl;
}
