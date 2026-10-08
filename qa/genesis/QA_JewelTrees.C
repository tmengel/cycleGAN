// QA / summary plots for JEWEL tree productions (Fun4All_JEWEL_MakeTrees.C output), in
// official sPHENIX plotting style. Overlays any number of samples.
//
// Reads, per sample:
//  - the flat towertree/jettree produced by Fun4All_JEWEL_MakeTrees.C (pure ROOT, no
//    Fun4All/coresoftware dependency)
//  - the truthtree produced by Fun4All_JEWEL_TruthQA.C: generator-level (pre-GEANT4) primary
//    particle kinematics, read straight from G4TruthInfo -- GEANT4 does not modify primary
//    particle 4-momenta, so this is equivalent to reading the HepMC files directly.
//
// Usage:
//   root -l -b -q 'QA_JewelTrees.C("<outdir>", "<sample spec>", "<energy label>")'
// with the sample spec a ';'-separated list of samples, each
//   name|legend label|ROOT color index|tree file glob|truthtree file
// The defaults reproduce the v1 (2025-08-22) vacuum/norecoil plots.
//
// Pages 12-15 (jet girth / N constituents) need the per-jet constituent branches
// (rjet_constit_*, tjet_constit_*); samples whose trees predate them are left out of those
// pages, and the pages are skipped entirely if no sample has them. Girth is
// g = sum_i (pT_i / pT_jet) * dR(i, jet), for jets with 20 < pT < 30 GeV and
// |eta| <= 1.1 - R = 0.7 -- the same definition and cuts as luisvale's
// plotGenesisObservables.C, so these pages are directly comparable to his.

#include "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.572/rootmacros/sPhenixStyle.C"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

struct Sample
{
  std::string name;
  std::string label;
  int color;
  TChain *towertree;
  TChain *jettree;
  TChain *truthtree;
  bool hasConstituents;
};

void DrawSPhenixLabel(double x, double y, const std::string &line2 = "")
{
  TLatex l;
  l.SetNDC();
  l.SetTextFont(62);
  l.SetTextSize(0.045);
  l.DrawLatex(x, y, "#it{#bf{sPHENIX}} Simulation");
  if (!line2.empty())
  {
    TLatex l2;
    l2.SetNDC();
    l2.SetTextFont(42);
    l2.SetTextSize(0.035);
    l2.DrawLatex(x, y - 0.055, line2.c_str());
  }
}

// Label for an individual pad in a multi-pad canvas (no legend possible on a colz plot).
void DrawPadLabel(const std::string &text)
{
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.045);
  l.DrawLatex(0.16, 0.925, text.c_str());  // in the pad's top margin, above the frame
}

std::vector<Sample> ParseSamples(const std::string &spec)
{
  std::vector<Sample> samples;
  TObjArray *entries = TString(spec).Tokenize(";");
  for (int i = 0; i < entries->GetEntries(); i++)
  {
    TObjArray *f = ((TObjString *) entries->At(i))->GetString().Tokenize("|");
    if (f->GetEntries() != 5)
    {
      std::cerr << "bad sample spec entry: " << ((TObjString *) entries->At(i))->GetString() << std::endl;
      continue;
    }
    auto field = [&](int k) { return std::string(((TObjString *) f->At(k))->GetString().Data()); };
    Sample s;
    s.name = field(0);
    s.label = field(1);
    s.color = std::stoi(field(2));
    s.towertree = new TChain("towertree");
    s.towertree->Add(field(3).c_str());
    s.jettree = new TChain("jettree");
    s.jettree->Add(field(3).c_str());
    s.truthtree = new TChain("truthtree");
    s.truthtree->Add(field(4).c_str());
    s.hasConstituents = s.jettree->GetBranch("rjet_constit_pt") != nullptr;
    samples.push_back(s);
  }
  return samples;
}

double JetGirth(float jpt, float jeta, float jphi, const std::vector<float> &cpt,
                const std::vector<float> &ceta, const std::vector<float> &cphi)
{
  if (jpt <= 0) return -1;
  double g = 0;
  for (size_t i = 0; i < cpt.size(); i++)
  {
    double deta = ceta[i] - jeta;
    double dphi = TVector2::Phi_mpi_pi(cphi[i] - jphi);
    g += cpt[i] * std::sqrt(deta * deta + dphi * dphi);
  }
  return g / jpt;
}

const std::string kDefaultSpec =
    "vacuum|JEWEL vacuum|862|/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/vacuum/trees/*.root|"
    "/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/qa_plots/truthqa_vacuum.root;"
    "norecoil|JEWEL no-recoil|633|/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/norecoil/trees/*.root|"
    "/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/qa_plots/truthqa_norecoil.root";

void QA_JewelTrees(const std::string &outdir = "/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/qa_plots",
                   const std::string &spec = kDefaultSpec,
                   const std::string &kEnergyLabel = "JEWEL p+p, #sqrt{s} = 200 GeV")
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kBird);

  gSystem->mkdir(outdir.c_str(), kTRUE);
  const std::string pdfname = outdir + "/JEWEL_QA_summary.pdf";

  std::vector<Sample> samples = ParseSamples(spec);
  if (samples.empty()) return;

  for (auto &s : samples)
  {
    std::cout << s.name << ": towertree=" << s.towertree->GetEntries()
              << " jettree=" << s.jettree->GetEntries()
              << " truthtree=" << s.truthtree->GetEntries()
              << " constituents=" << (s.hasConstituents ? "yes" : "no") << std::endl;
  }

  // canvas grid for per-sample 2D pages
  const int ncol = samples.size() <= 2 ? samples.size() : 2;
  const int nrow = (samples.size() + ncol - 1) / ncol;

  TCanvas *c = new TCanvas("c", "c", 900, 700);
  c->Print((pdfname + "[").c_str());

  // Draws prepared histograms as an overlay. Legend and sPHENIX label sit in the headroom
  // above the tallest bin (x50 on log scale, x1.4 on linear), which is empty at every x, so
  // they can't land on the data whatever the spectrum shape.
  auto drawOverlay = [&](const std::string &pngname, std::vector<TH1D *> &hists,
                         const std::vector<std::string> &labels, bool logy, const std::string &extraLabel,
                         const std::string &cutLabel = "")
  {
    c->Clear();
    c->SetLogy(logy);
    // right-aligned against the frame edge, clear of the top-left label block
    TLegend *leg = new TLegend(0.55, 0.91 - 0.048 * hists.size(), 0.94, 0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.028);
    leg->SetMargin(0.18);
    double ymax = 0;
    double ymin = 1e30;
    for (size_t i = 0; i < hists.size(); i++)
    {
      ymax = std::max(ymax, hists[i]->GetMaximum());
      for (int b = 1; b <= hists[i]->GetNbinsX(); b++)
      {
        if (hists[i]->GetBinContent(b) > 0) ymin = std::min(ymin, hists[i]->GetBinContent(b));
      }
      hists[i]->GetYaxis()->SetTitleOffset(1.6);
      leg->AddEntry(hists[i], labels[i].c_str(), "l");
    }
    // linear plots with a label get extra headroom (x1.9) so the label block fits above the data
    hists[0]->SetMaximum(logy ? ymax * 50 : ymax * (extraLabel.empty() ? 1.4 : 1.9));
    hists[0]->SetMinimum(logy ? std::max(ymin * 0.5, ymax * 1e-7) : 0);
    hists[0]->Draw("hist");
    for (size_t i = 1; i < hists.size(); i++) hists[i]->Draw("hist same");
    if (!extraLabel.empty()) DrawSPhenixLabel(0.19, 0.86, extraLabel);
    if (!cutLabel.empty())
    {
      TLatex l3;
      l3.SetNDC();
      l3.SetTextFont(42);
      l3.SetTextSize(0.035);
      l3.DrawLatex(0.19, 0.86 - 0.11, cutLabel.c_str());
    }
    leg->Draw();
    c->Print(pdfname.c_str());
    c->SaveAs((outdir + "/" + pngname + ".png").c_str());
  };

  // density=true (default) plots a proper (1/N_evt) dN/dX differential -- also divides out
  // the bin width, so the y-axis is a true per-unit-X density and stays correct regardless of
  // binning. Set false for a histogram whose x-axis is already a discrete per-event count
  // (e.g. multiplicity) where "per unit X" isn't a meaningful density and the y-axis is just
  // a fraction of events.
  auto overlay1D = [&](const std::string &pngname, const char *drawexpr, const char *cut,
                        const char *title, int nbin, double lo, double hi, bool logy,
                        TChain *Sample::*tree, const std::string &extraLabel = "", bool density = true)
  {
    std::vector<TH1D *> hists;
    std::vector<std::string> labels;
    for (auto &s : samples)
    {
      TChain *t = s.*tree;
      TH1D *h = new TH1D(Form("h_%s_%s", pngname.c_str(), s.name.c_str()), title, nbin, lo, hi);
      t->Draw(Form("%s>>%s", drawexpr, h->GetName()), cut, "goff");
      if (h->Integral() > 0) h->Scale(1. / t->GetEntries());
      if (density) h->Scale(1.0 / h->GetBinWidth(1));
      h->SetLineColor(s.color);
      h->SetLineWidth(2);
      hists.push_back(h);
      labels.push_back(s.label);
    }
    drawOverlay(pngname, hists, labels, logy, extraLabel);
  };

  // Multi-pad overlay: one pad per branch, all samples overlaid in each pad.
  auto overlayPads = [&](const std::string &pngname, const std::vector<const char *> &branches,
                         const std::vector<const char *> &titles, int nbin, double lo, double hi,
                         TChain *Sample::*tree, bool legendTop)
  {
    c->Clear();
    c->SetLogy(0);
    c->Divide(branches.size(), 1);
    for (size_t pad = 0; pad < branches.size(); pad++)
    {
      c->cd(pad + 1);
      gPad->SetLeftMargin(0.19);
      TLegend *leg = legendTop ? new TLegend(0.22, 0.90 - 0.06 * samples.size(), 0.62, 0.90)
                               : new TLegend(0.22, 0.18, 0.62, 0.18 + 0.06 * samples.size());
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.035);
      std::vector<TH1D *> hists;
      for (auto &s : samples)
      {
        TChain *t = s.*tree;
        TH1D *h = new TH1D(Form("h_%s_%zu_%s", pngname.c_str(), pad, s.name.c_str()), titles[pad], nbin, lo, hi);
        t->Draw(Form("%s>>%s", branches[pad], h->GetName()), "", "goff");
        if (t->GetEntries() > 0) h->Scale(1. / t->GetEntries());
        h->Scale(1.0 / h->GetBinWidth(1));
        h->SetLineColor(s.color);
        h->SetLineWidth(2);
        h->SetMinimum(0);
        h->GetYaxis()->SetTitleOffset(2.1);
        hists.push_back(h);
        leg->AddEntry(h, s.label.c_str(), "l");
      }
      double ymax = 0;
      for (auto h : hists) ymax = std::max(ymax, h->GetMaximum());
      // legend band on top: headroom grows with the number of legend rows (x1.6 for 2, x2.1 for 4)
      hists[0]->SetMaximum(ymax * (legendTop ? std::max(1.6, 1.1 + 0.25 * samples.size()) : 1.4));
      hists[0]->Draw("hist");
      for (size_t i = 1; i < hists.size(); i++) hists[i]->Draw("hist same");
      leg->Draw();
    }
    c->Print(pdfname.c_str());
    c->SaveAs((outdir + "/" + pngname + ".png").c_str());
  };

  //======================================================================
  // Page 1: truth-level (pre-GEANT4) primary particle pT spectrum
  //======================================================================
  overlay1D("page1_truthlevel_particle_pt", "truth_pt", "1",
            ";generator-level particle p_{T} [GeV];(1/N_{evt}) dN/dp_{T} [GeV^{-1}]",
            40, 0, 8, true, &Sample::truthtree, kEnergyLabel);

  //======================================================================
  // Page 2: truth-level particle eta / phi (2-pad)
  //======================================================================
  {
    c->Clear();
    c->SetLogy(0);
    c->Divide(2, 1);
    for (int pad = 0; pad < 2; pad++)
    {
      c->cd(pad + 1);
      gPad->SetLeftMargin(0.19);
      TLegend *leg = new TLegend(0.22, 0.90 - 0.06 * samples.size(), 0.62, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.035);
      std::vector<TH1D *> hists;
      const char *branch = pad == 0 ? "truth_eta" : "truth_phi";
      const char *title = pad == 0 ? ";generator-level particle #eta;(1/N_{evt}) dN/d#eta"
                                   : ";generator-level particle #phi [rad];(1/N_{evt}) dN/d#phi";
      double lo = pad == 0 ? -5 : -TMath::Pi();
      double hi = pad == 0 ? 5 : TMath::Pi();
      for (auto &s : samples)
      {
        TH1D *h = new TH1D(Form("h_truth_%s_%s", branch, s.name.c_str()), title, 40, lo, hi);
        s.truthtree->Draw(Form("%s>>%s", branch, h->GetName()), "", "goff");
        if (s.truthtree->GetEntries() > 0) h->Scale(1. / s.truthtree->GetEntries());
        h->Scale(1.0 / h->GetBinWidth(1));
        h->SetLineColor(s.color);
        h->SetLineWidth(2);
        h->SetMinimum(0);
        h->GetYaxis()->SetTitleOffset(2.1);
        hists.push_back(h);
        leg->AddEntry(h, s.label.c_str(), "l");
      }
      double ymax = 0;
      for (auto h : hists) ymax = std::max(ymax, h->GetMaximum());
      hists[0]->SetMaximum(ymax * std::max(1.6, 1.1 + 0.25 * samples.size()));
      hists[0]->Draw("hist");
      for (size_t i = 1; i < hists.size(); i++) hists[i]->Draw("hist same");
      leg->Draw();
    }
    c->Print(pdfname.c_str());
    c->SaveAs((outdir + "/page2_truthlevel_particle_eta_phi.png").c_str());
  }

  //======================================================================
  // Page 3: truth-level particle multiplicity per event
  //======================================================================
  overlay1D("page3_truthlevel_multiplicity", "Length$(truth_pt)", "1",
            ";generator-level particles / event;Fraction of events", 40, 0, 200, false,
            &Sample::truthtree, "", false);

  //======================================================================
  // Page 4: leading truth jet pT (post-GEANT4 reco chain)
  //======================================================================
  // axis starts at the generator jet filter threshold (20 GeV): below it the spectrum is empty
  overlay1D("page4_leading_truth_jet_pt", "Max$(tjet_pt)", "Length$(tjet_pt)>0",
            ";leading truth jet p_{T} [GeV];(1/N_{evt}) dN/dp_{T} [GeV^{-1}]", 32, 18, 66, true,
            &Sample::jettree, kEnergyLabel);

  //======================================================================
  // Page 5: leading reco jet pT
  //======================================================================
  // axis starts at the tree's own jet cut (10 GeV)
  overlay1D("page5_leading_reco_jet_pt", "Max$(rjet_pt)", "Length$(rjet_pt)>0",
            ";leading reco jet p_{T} [GeV];(1/N_{evt}) dN/dp_{T} [GeV^{-1}]", 28, 10, 66, true,
            &Sample::jettree, kEnergyLabel);

  //======================================================================
  // Pages 6-7: jet eta, jet phi (truth, reco) -- 2-pad
  //======================================================================
  overlayPads("page6_jet_eta", {"tjet_eta", "rjet_eta"},
              {";truth jet #eta;(1/N_{evt}) dN/d#eta", ";reco jet #eta;(1/N_{evt}) dN/d#eta"},
              24, -1.2, 1.2, &Sample::jettree, true);
  overlayPads("page7_jet_phi", {"tjet_phi", "rjet_phi"},
              {";truth jet #phi [rad];(1/N_{evt}) dN/d#phi", ";reco jet #phi [rad];(1/N_{evt}) dN/d#phi"},
              32, -TMath::Pi(), TMath::Pi(), &Sample::jettree, true);

  //======================================================================
  // Page 8: jet multiplicity per event (truth, reco) -- 2-pad
  //======================================================================
  {
    c->Clear();
    c->Divide(2, 1);
    for (int pad = 0; pad < 2; pad++)
    {
      c->cd(pad + 1);
      gPad->SetLeftMargin(0.16);
      gPad->SetLogy();
      TLegend *leg = new TLegend(0.50, 0.90 - 0.06 * samples.size(), 0.90, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.03);
      std::vector<TH1D *> hists;
      const char *branch = pad == 0 ? "Length$(tjet_pt)" : "Length$(rjet_pt)";
      const char *title = pad == 0 ? ";N truth jets / event;Fraction of events" : ";N reco jets / event;Fraction of events";
      for (auto &s : samples)
      {
        TH1D *h = new TH1D(Form("h_njet_%d_%s", pad, s.name.c_str()), title, 10, 0, 10);
        s.jettree->Draw(Form("%s>>%s", branch, h->GetName()), "", "goff");
        h->Scale(1. / s.jettree->GetEntries());
        h->SetLineColor(s.color);
        h->SetLineWidth(2);
        h->GetYaxis()->SetTitleOffset(1.8);
        hists.push_back(h);
        leg->AddEntry(h, s.label.c_str(), "l");
      }
      hists[0]->SetMaximum(50);
      hists[0]->SetMinimum(1e-5);
      hists[0]->Draw("hist");
      for (size_t i = 1; i < hists.size(); i++) hists[i]->Draw("hist same");
      leg->Draw();
    }
    c->Print(pdfname.c_str());
    c->SaveAs((outdir + "/page8_jet_multiplicity.png").c_str());
  }

  //======================================================================
  // Page 9: 2D truth vs reco leading jet pT correlation, one pad per sample
  //======================================================================
  {
    // own canvas sized to the pad grid: a shared 900x700 canvas clips the z-axis titles
    TCanvas *c2d = new TCanvas("c_page9", "c_page9", 700 * ncol, 600 * nrow);
    c2d->Divide(ncol, nrow);
    for (size_t i = 0; i < samples.size(); i++)
    {
      c2d->cd(i + 1);
      gPad->SetTopMargin(0.10);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.25);
      gPad->SetLogz();
      auto &s = samples[i];
      TH2D *h2 = new TH2D(Form("h2_truthvsreco_%s", s.name.c_str()),
                          ";leading truth jet p_{T} [GeV];leading reco jet p_{T} [GeV]",
                          30, 0, 60, 30, 0, 60);
      s.jettree->Draw(Form("Max$(rjet_pt):Max$(tjet_pt)>>%s", h2->GetName()),
                      "Length$(tjet_pt)>0 && Length$(rjet_pt)>0", "goff");
      h2->Scale(1. / s.jettree->GetEntries());
      h2->GetZaxis()->SetTitle("Fraction of events");
      h2->GetZaxis()->SetTitleOffset(1.75);
      h2->Draw("colz");
      TLine line;
      line.SetLineStyle(2);
      line.DrawLine(0, 0, 60, 60);
      DrawPadLabel(s.label);
    }
    c2d->Print(pdfname.c_str());
    c2d->SaveAs((outdir + "/page9_truth_vs_reco_pt.png").c_str());
  }

  //======================================================================
  // Page 10: average calo image (eta-phi energy map), one pad per sample
  //======================================================================
  {
    // own canvas sized to the pad grid: a shared 900x700 canvas clips the z-axis titles
    TCanvas *c2d = new TCanvas("c_page10", "c_page10", 700 * ncol, 600 * nrow);
    c2d->Divide(ncol, nrow);
    double pi = TMath::Pi();
    for (size_t i = 0; i < samples.size(); i++)
    {
      c2d->cd(i + 1);
      gPad->SetTopMargin(0.10);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.25);
      auto &s = samples[i];
      TH2D *h2avg = new TH2D(Form("h2_avgimg_%s", s.name.c_str()), ";#eta;#phi [rad]",
                             24, -1.1, 1.1, 64, 0, 2 * pi);

      std::vector<float> *em = nullptr, *ih = nullptr, *oh = nullptr;
      std::vector<int> *etabin = nullptr, *phibin = nullptr;
      s.towertree->SetBranchStatus("*", 0);
      for (auto b : {"emcalenergy", "ihcalenergy", "ohcalenergy", "etabin", "phibin"}) s.towertree->SetBranchStatus(b, 1);
      s.towertree->SetBranchAddress("emcalenergy", &em);
      s.towertree->SetBranchAddress("ihcalenergy", &ih);
      s.towertree->SetBranchAddress("ohcalenergy", &oh);
      s.towertree->SetBranchAddress("etabin", &etabin);
      s.towertree->SetBranchAddress("phibin", &phibin);

      Long64_t n = s.towertree->GetEntries();
      for (Long64_t e = 0; e < n; e++)
      {
        s.towertree->GetEntry(e);
        int ncal = em->size();
        for (int k = 0; k < ncal; k++)
        {
          double tot = std::max(0.f, (*em)[k]) + std::max(0.f, (*ih)[k]) + std::max(0.f, (*oh)[k]);
          h2avg->SetBinContent((*etabin)[k] + 1, (*phibin)[k] + 1,
                               h2avg->GetBinContent((*etabin)[k] + 1, (*phibin)[k] + 1) + tot);
        }
      }
      if (n > 0) h2avg->Scale(1. / n);
      h2avg->GetZaxis()->SetTitle("#LT E #GT [GeV]");
      h2avg->GetZaxis()->SetTitleOffset(1.75);
      h2avg->Draw("colz");
      DrawPadLabel(s.label);
      s.towertree->ResetBranchAddresses();
      s.towertree->SetBranchStatus("*", 1);
    }
    c2d->Print(pdfname.c_str());
    c2d->SaveAs((outdir + "/page10_average_calo_image.png").c_str());
  }

  //======================================================================
  // Page 11: per-event total calo energy
  //======================================================================
  overlay1D("page11_total_calo_energy",
            "Sum$(max(emcalenergy,0)+max(ihcalenergy,0)+max(ohcalenergy,0))", "1",
            ";#Sigma E_{tower} per event [GeV];(1/N_{evt}) dN/dE [GeV^{-1}]", 40, 0, 200, true,
            &Sample::towertree, kEnergyLabel);

  //======================================================================
  // Pages 12-15: jet substructure from constituents (20 < pT < 30 GeV, |eta| <= 0.7)
  //======================================================================
  const float kR = 0.4;
  const float kPtLo = 20, kPtHi = 30;
  std::vector<Sample *> withConst;
  for (auto &s : samples)
  {
    if (s.hasConstituents) withConst.push_back(&s);
  }
  if (withConst.empty())
  {
    std::cout << "no sample has constituent branches -- skipping girth / N-constituent pages" << std::endl;
  }
  else
  {
    std::vector<TH1D *> hGirthR, hGirthT, hNcR, hNcT;
    std::vector<std::string> labels;
    for (Sample *s : withConst)
    {
      TH1D *gR = new TH1D(Form("hGirthReco_%s", s->name.c_str()), ";reco jet girth g;(1/N_{jet}) dN/dg", 40, 0, kR);
      TH1D *gT = new TH1D(Form("hGirthTruth_%s", s->name.c_str()), ";truth jet girth g;(1/N_{jet}) dN/dg", 40, 0, kR);
      TH1D *nR = new TH1D(Form("hNConstReco_%s", s->name.c_str()), ";reco jet N constituents (towers);Fraction of jets", 60, -0.5, 119.5);
      TH1D *nT = new TH1D(Form("hNConstTruth_%s", s->name.c_str()), ";truth jet N constituents (particles);Fraction of jets", 40, -0.5, 39.5);

      std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr;
      std::vector<std::vector<float>> *cpt = nullptr, *ceta = nullptr, *cphi = nullptr;
      TChain *t = s->jettree;
      for (int pass = 0; pass < 2; pass++)
      {
        const std::string pre = pass == 0 ? "rjet" : "tjet";
        TH1D *hg = pass == 0 ? gR : gT;
        TH1D *hn = pass == 0 ? nR : nT;
        t->SetBranchStatus("*", 0);
        for (auto b : {"_pt", "_eta", "_phi", "_constit_pt", "_constit_eta", "_constit_phi"}) t->SetBranchStatus((pre + b).c_str(), 1);
        t->SetBranchAddress((pre + "_pt").c_str(), &pt);
        t->SetBranchAddress((pre + "_eta").c_str(), &eta);
        t->SetBranchAddress((pre + "_phi").c_str(), &phi);
        t->SetBranchAddress((pre + "_constit_pt").c_str(), &cpt);
        t->SetBranchAddress((pre + "_constit_eta").c_str(), &ceta);
        t->SetBranchAddress((pre + "_constit_phi").c_str(), &cphi);
        Long64_t n = t->GetEntries();
        for (Long64_t e = 0; e < n; e++)
        {
          t->GetEntry(e);
          for (size_t j = 0; j < pt->size(); j++)
          {
            if ((*pt)[j] < kPtLo || (*pt)[j] > kPtHi || std::fabs((*eta)[j]) > 1.1 - kR) continue;
            hn->Fill((*cpt)[j].size());
            double g = JetGirth((*pt)[j], (*eta)[j], (*phi)[j], (*cpt)[j], (*ceta)[j], (*cphi)[j]);
            if (g >= 0) hg->Fill(g);
          }
        }
        t->ResetBranchAddresses();
        t->SetBranchStatus("*", 1);
      }
      for (TH1D *h : {gR, gT, nR, nT})
      {
        if (h->GetEntries() > 0) h->Scale(1. / h->GetEntries());
        h->SetLineColor(s->color);
        h->SetLineWidth(2);
      }
      gR->Scale(1. / gR->GetBinWidth(1));
      gT->Scale(1. / gT->GetBinWidth(1));
      hGirthR.push_back(gR);
      hGirthT.push_back(gT);
      hNcR.push_back(nR);
      hNcT.push_back(nT);
      labels.push_back(s->label);
    }
    const std::string cutLabel = "anti-k_{T} R=0.4, 20 < p_{T}^{jet} < 30 GeV, |#eta^{jet}| < 0.7";
    drawOverlay("page12_girth_truth", hGirthT, labels, false, kEnergyLabel, cutLabel);
    drawOverlay("page13_girth_reco", hGirthR, labels, false, kEnergyLabel, cutLabel);
    drawOverlay("page14_nconstituents_truth", hNcT, labels, false, kEnergyLabel, cutLabel);
    drawOverlay("page15_nconstituents_reco", hNcR, labels, false, kEnergyLabel, cutLabel);
  }

  c->Print((pdfname + "]").c_str());
  std::cout << "Wrote " << pdfname << std::endl;
}
