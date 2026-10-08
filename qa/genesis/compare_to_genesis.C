// Compares luisvale's genesisObservables_R4.root (norecoil/vacuum, R=0.4, DijetTreeMaker-based
// production from plotGenesisObservables.C) against the same observables computed from a
// tmengel production (Fun4All_JEWEL_MakeTrees.C trees + Fun4All_JEWEL_TruthQA.C truthtree).
//
// Usage:
//   root -l -b -q 'compare_to_genesis.C+("<tree base dir>", "<truthqa dir>", "<outdir>",
//                                        "<genesis file>", "<production legend label>")'
// reading <tree base dir>/<sim>/trees/*.root and <truthqa dir>/truthqa_<sim>.root.
// Defaults reproduce the v1 comparison paths.
//
// Every observable is booked with genesis's exact binning and filled with genesis's cuts:
//  - jet eta/phi, leading-jet pT, N jets, girth, N constituents, truth-jet constituent
//    kinematics: |eta_jet| <= 1.1 - R = 0.7
//  - girth / N constituents: additionally 20 <= pT_jet <= 30 GeV;
//    girth g = sum_i (pT_i / pT_jet) * dR(i, jet)  (plotGenesisObservables.C's jetGirthG)
//  - hGenPart*: constituents of fiducial truth jets (genesis's "generator-level" proxy)
//  - hHepMCPart*: all generator-level particles (our truthtree = G4TruthInfo primaries,
//    which GEANT4 leaves untouched, vs genesis's status==1 HepMC particles)
//
// Caveats (not fixable here):
//  - our jettree jets (rjet_*, tjet_*) are cut at pT > 10 GeV at write time
//    (JewelTreeWriter m_jet_pt_min); a dashed line at 10 GeV marks this on pT plots.
//  - reco jet constituents are calorimeter towers in both productions, but tower inputs
//    (retowered EMCal + IHCal + OHCal) and noise treatment may differ from DijetTreeMaker's,
//    so reco N-constituents in particular is a chain comparison, not a closure test.
//  - trees without constituent branches (v1) skip the girth / N-constituent / hGenPart plots.
//  - genesisObservables_R4.root is dated 2025-08-25, i.e. made from luisvale's earlier
//    10k-event merge. His vacuum JEWEL parameters did not change after that, but his norecoil
//    medium did (TI 0.260 -> 0.375 GeV, medium_auau_c0.params.dat, 2025-09-01). A newer
//    production's norecoil comparison therefore includes that physics change.

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
#include <TPad.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  const float R = 0.4;
  const float ETAFID = 1.1 - R;  // 0.7
  const float PTLO = 20, PTHI = 30;

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
}  // namespace

struct ProdHists
{
  TH1D *hJetPtReco = nullptr, *hJetPtTruth = nullptr;
  TH1D *hJetEta = nullptr, *hJetPhi = nullptr;
  TH1D *hLeadJetPtReco = nullptr, *hLeadJetPtTruth = nullptr;
  TH1D *hNJetsReco = nullptr, *hNJetsTruth = nullptr;
  TH2D *hLeadTruthVsRecoPt = nullptr;
  TH1D *hHepMCPartPt = nullptr, *hHepMCPartEta = nullptr, *hHepMCPartPhi = nullptr;
  bool hasConstituents = false;
  TH1D *hNConstitReco = nullptr, *hNConstitTruth = nullptr;
  TH1D *hGirthReco = nullptr, *hGirthTruth = nullptr;
  TProfile *pGirthVsPt = nullptr;
  TH1D *hGenPartPt = nullptr, *hGenPartEta = nullptr, *hGenPartPhi = nullptr;
};

// Builds production histograms for one sim, binning matched to genesis's histogram of the
// same name, with the same cuts.
ProdHists buildProdHists(const std::string &treebase, const std::string &truthdir, const char *sim)
{
  ProdHists h;

  TChain jettree("jettree");
  jettree.Add(Form("%s/%s/trees/*.root", treebase.c_str(), sim));
  TChain truthtree("truthtree");
  truthtree.Add(Form("%s/truthqa_%s.root", truthdir.c_str(), sim));
  h.hasConstituents = jettree.GetBranch("rjet_constit_pt") != nullptr;

  std::vector<float> *rjet_pt = nullptr, *rjet_eta = nullptr, *rjet_phi = nullptr;
  std::vector<float> *tjet_pt = nullptr, *tjet_eta = nullptr, *tjet_phi = nullptr;
  std::vector<std::vector<float>> *rc_pt = nullptr, *rc_eta = nullptr, *rc_phi = nullptr;
  std::vector<std::vector<float>> *tc_pt = nullptr, *tc_eta = nullptr, *tc_phi = nullptr;
  jettree.SetBranchAddress("rjet_pt", &rjet_pt);
  jettree.SetBranchAddress("rjet_eta", &rjet_eta);
  jettree.SetBranchAddress("rjet_phi", &rjet_phi);
  jettree.SetBranchAddress("tjet_pt", &tjet_pt);
  jettree.SetBranchAddress("tjet_eta", &tjet_eta);
  jettree.SetBranchAddress("tjet_phi", &tjet_phi);
  if (h.hasConstituents)
  {
    jettree.SetBranchAddress("rjet_constit_pt", &rc_pt);
    jettree.SetBranchAddress("rjet_constit_eta", &rc_eta);
    jettree.SetBranchAddress("rjet_constit_phi", &rc_phi);
    jettree.SetBranchAddress("tjet_constit_pt", &tc_pt);
    jettree.SetBranchAddress("tjet_constit_eta", &tc_eta);
    jettree.SetBranchAddress("tjet_constit_phi", &tc_phi);
  }

  h.hJetPtReco = new TH1D(Form("p_hJetPtReco_%s", sim), ";reco jet p_{T} [GeV]", 100, 0, 100);
  h.hJetPtTruth = new TH1D(Form("p_hJetPtTruth_%s", sim), ";truth jet p_{T} [GeV]", 100, 0, 100);
  h.hJetEta = new TH1D(Form("p_hJetEta_%s", sim), ";reco jet #eta", 50, -1.2, 1.2);
  h.hJetPhi = new TH1D(Form("p_hJetPhi_%s", sim), ";reco jet #phi [rad]", 50, -TMath::Pi(), TMath::Pi());
  h.hLeadJetPtReco = new TH1D(Form("p_hLeadJetPtReco_%s", sim), ";leading reco jet p_{T} [GeV]", 100, 0, 100);
  h.hLeadJetPtTruth = new TH1D(Form("p_hLeadJetPtTruth_%s", sim), ";leading truth jet p_{T} [GeV]", 100, 0, 100);
  h.hNJetsReco = new TH1D(Form("p_hNJetsReco_%s", sim), ";N reco jets / event", 15, -0.5, 14.5);
  h.hNJetsTruth = new TH1D(Form("p_hNJetsTruth_%s", sim), ";N truth jets / event", 15, -0.5, 14.5);
  // genesis fills (x = reco, y = truth) -- matched here
  h.hLeadTruthVsRecoPt = new TH2D(Form("p_hLeadTruthVsRecoPt_%s", sim), "", 60, 0, 60, 60, 0, 60);
  h.hNConstitReco = new TH1D(Form("p_hNConstitReco_%s", sim), ";reco jet N constituents", 125, -0.5, 249.5);
  h.hNConstitTruth = new TH1D(Form("p_hNConstitTruth_%s", sim), ";truth jet N constituents", 40, -0.5, 39.5);
  h.hGirthReco = new TH1D(Form("p_hGirthReco_%s", sim), ";reco jet girth g", 100, 0, R);
  h.hGirthTruth = new TH1D(Form("p_hGirthTruth_%s", sim), ";truth jet girth g", 100, 0, R);
  h.pGirthVsPt = new TProfile(Form("p_pGirthVsPt_%s", sim), ";reco jet p_{T} [GeV];#LT g #GT", 50, 0, 100);
  h.hGenPartPt = new TH1D(Form("p_hGenPartPt_%s", sim), ";truth-jet constituent p_{T} [GeV]", 100, 0, 30);
  h.hGenPartEta = new TH1D(Form("p_hGenPartEta_%s", sim), ";truth-jet constituent #eta", 50, -1.2, 1.2);
  h.hGenPartPhi = new TH1D(Form("p_hGenPartPhi_%s", sim), ";truth-jet constituent #phi [rad]", 50, -TMath::Pi(), TMath::Pi());

  Long64_t nEntries = jettree.GetEntries();
  for (Long64_t i = 0; i < nEntries; i++)
  {
    jettree.GetEntry(i);
    float leadReco = -1, leadTruth = -1;
    int nRecoFid = 0, nTruthFid = 0;
    for (size_t j = 0; j < rjet_pt->size(); j++)
    {
      float pt = rjet_pt->at(j), eta = rjet_eta->at(j), phi = rjet_phi->at(j);
      h.hJetPtReco->Fill(pt);
      if (std::fabs(eta) > ETAFID) continue;
      h.hJetEta->Fill(eta);
      h.hJetPhi->Fill(phi);
      nRecoFid++;
      if (pt > leadReco) leadReco = pt;
      if (!h.hasConstituents) continue;
      double g = JetGirth(pt, eta, phi, rc_pt->at(j), rc_eta->at(j), rc_phi->at(j));
      if (g >= 0) h.pGirthVsPt->Fill(pt, g);
      if (pt >= PTLO && pt <= PTHI)
      {
        h.hNConstitReco->Fill(rc_pt->at(j).size());
        if (g >= 0) h.hGirthReco->Fill(g);
      }
    }
    h.hNJetsReco->Fill(nRecoFid);
    if (leadReco >= 0) h.hLeadJetPtReco->Fill(leadReco);

    for (size_t j = 0; j < tjet_pt->size(); j++)
    {
      float pt = tjet_pt->at(j), eta = tjet_eta->at(j), phi = tjet_phi->at(j);
      h.hJetPtTruth->Fill(pt);
      if (std::fabs(eta) > ETAFID) continue;
      nTruthFid++;
      if (pt > leadTruth) leadTruth = pt;
      if (!h.hasConstituents) continue;
      const std::vector<float> &cpt = tc_pt->at(j);
      for (size_t ic = 0; ic < cpt.size(); ic++)
      {
        h.hGenPartPt->Fill(cpt[ic]);
        h.hGenPartEta->Fill(tc_eta->at(j)[ic]);
        h.hGenPartPhi->Fill(tc_phi->at(j)[ic]);
      }
      if (pt >= PTLO && pt <= PTHI)
      {
        h.hNConstitTruth->Fill(cpt.size());
        double g = JetGirth(pt, eta, phi, cpt, tc_eta->at(j), tc_phi->at(j));
        if (g >= 0) h.hGirthTruth->Fill(g);
      }
    }
    h.hNJetsTruth->Fill(nTruthFid);
    if (leadTruth >= 0) h.hLeadJetPtTruth->Fill(leadTruth);
    if (leadReco >= 0 && leadTruth >= 0) h.hLeadTruthVsRecoPt->Fill(leadReco, leadTruth);
  }

  h.hHepMCPartPt = new TH1D(Form("p_hHepMCPartPt_%s", sim), ";generator-level particle p_{T} [GeV]", 150, 0, 30);
  h.hHepMCPartEta = new TH1D(Form("p_hHepMCPartEta_%s", sim), ";generator-level particle #eta", 120, -10, 10);
  h.hHepMCPartPhi = new TH1D(Form("p_hHepMCPartPhi_%s", sim), ";generator-level particle #phi [rad]", 50, -TMath::Pi(), TMath::Pi());
  // TTree::Draw(">>name") finds its target histogram by name in gDirectory. With
  // TH1::AddDirectory(kFALSE) these aren't registered anywhere, and looping over the jettree
  // chain has left gDirectory pointing at its last file -- so cd to gROOT and attach them there
  gROOT->cd();
  for (TH1D *hh : {h.hHepMCPartPt, h.hHepMCPartEta, h.hHepMCPartPhi}) hh->SetDirectory(gROOT);
  truthtree.Draw(Form("truth_pt>>p_hHepMCPartPt_%s", sim), "", "goff");
  truthtree.Draw(Form("truth_eta>>p_hHepMCPartEta_%s", sim), "", "goff");
  truthtree.Draw(Form("truth_phi>>p_hHepMCPartPhi_%s", sim), "", "goff");
  for (TH1D *hh : {h.hHepMCPartPt, h.hHepMCPartEta, h.hHepMCPartPhi}) hh->SetDirectory(nullptr);

  std::cout << sim << " (production): " << nEntries << " events, "
            << h.hJetPtReco->GetEntries() << " reco jets, "
            << h.hJetPtTruth->GetEntries() << " truth jets, "
            << truthtree.GetEntries() << " truthtree events, constituents="
            << (h.hasConstituents ? "yes" : "no") << std::endl;

  return h;
}

// Shape comparison (both normalized to unit area) with a ratio panel (production / genesis).
void drawOverlay(TH1 *hGen, TH1 *hProd, const std::string &prodLabel, const std::string &label,
                 const std::string &outpng, bool logy, bool markTenGeV = false)
{
  if (!hGen || !hProd)
  {
    std::cerr << "missing histogram for " << outpng << std::endl;
    return;
  }
  TCanvas *c = new TCanvas("c_ov", "c_ov", 800, 800);
  TPad *top = new TPad("top", "top", 0, 0.30, 1, 1);
  TPad *bot = new TPad("bot", "bot", 0, 0, 1, 0.30);
  top->SetBottomMargin(0.02);
  top->SetLeftMargin(0.16);
  bot->SetTopMargin(0.03);
  bot->SetBottomMargin(0.35);
  bot->SetLeftMargin(0.16);
  top->Draw();
  bot->Draw();

  TH1 *g = (TH1 *) hGen->Clone("g_clone");
  TH1 *p = (TH1 *) hProd->Clone("p_clone");
  if (g->GetSumw2N() == 0) g->Sumw2();
  if (p->GetSumw2N() == 0) p->Sumw2();
  // pT plots: our jets are cut at 10 GeV at write time, genesis's are not, so normalize both
  // over pT > 10 GeV only -- otherwise genesis's low-pT jets would distort the shape comparison
  int b0 = markTenGeV ? g->GetXaxis()->FindBin(10.0 + 1e-3) : 1;
  int b1 = g->GetNbinsX();
  if (g->Integral(b0, b1) > 0) g->Scale(1.0 / g->Integral(b0, b1));
  if (p->Integral(b0, b1) > 0) p->Scale(1.0 / p->Integral(b0, b1));
  g->SetLineColor(kAzure + 2);
  g->SetMarkerColor(kAzure + 2);
  g->SetLineWidth(2);
  p->SetLineColor(kRed + 1);
  p->SetMarkerColor(kRed + 1);
  p->SetMarkerStyle(20);
  p->SetMarkerSize(0.7);

  top->cd();
  top->SetLogy(logy);
  g->GetXaxis()->SetLabelSize(0);
  g->GetXaxis()->SetTitleSize(0);
  g->GetYaxis()->SetTitle(markTenGeV ? "Normalized entries (p_{T} > 10 GeV)" : "Normalized entries");
  g->GetYaxis()->SetTitleOffset(1.5);
  double ymax = std::max(g->GetMaximum(), p->GetMaximum());
  double ymin = 1e30;
  for (int b = 1; b <= g->GetNbinsX(); b++)
  {
    if (g->GetBinContent(b) > 0) ymin = std::min(ymin, g->GetBinContent(b));
    if (p->GetBinContent(b) > 0) ymin = std::min(ymin, p->GetBinContent(b));
  }
  // headroom above the tallest bin keeps the label block and legend off the data
  g->SetMaximum(logy ? ymax * 3000 : ymax * 2.2);
  g->SetMinimum(logy ? std::max(ymin * 0.5, ymax * 1e-7) : 0);
  g->Draw("hist");
  p->Draw("e1 same");
  if (markTenGeV)
  {
    TLine *l = new TLine(10, g->GetMinimum(), 10, logy ? ymax * 3 : ymax * 1.1);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw();
  }
  // label lines (newline-separated) top-left, legend top-right beside the sPHENIX line only
  TLatex lt;
  lt.SetNDC();
  lt.SetTextFont(62);
  lt.SetTextSize(0.045);
  lt.DrawLatex(0.20, 0.86, "#it{#bf{sPHENIX}} Simulation");
  lt.SetTextFont(42);
  lt.SetTextSize(0.034);
  TObjArray *lines = TString(label).Tokenize("\n");
  for (int i = 0; i < lines->GetEntries(); i++)
  {
    lt.DrawLatex(0.20, 0.80 - 0.05 * i, ((TObjString *) lines->At(i))->GetString().Data());
  }
  TLegend *leg = new TLegend(0.60, 0.80, 0.93, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.034);
  leg->AddEntry(g, "genesis (luisvale)", "l");
  leg->AddEntry(p, prodLabel.c_str(), "pe");
  leg->Draw();

  bot->cd();
  TH1 *r = (TH1 *) p->Clone("r_clone");
  r->Divide(g);
  r->SetMinimum(0.0);
  r->SetMaximum(2.0);
  r->GetYaxis()->SetTitle("prod / genesis");
  r->GetYaxis()->SetNdivisions(505);
  r->GetYaxis()->SetTitleSize(0.10);
  r->GetYaxis()->SetTitleOffset(0.7);
  r->GetYaxis()->SetLabelSize(0.09);
  r->GetXaxis()->SetTitle(hProd->GetXaxis()->GetTitle());
  r->GetXaxis()->SetTitleSize(0.12);
  r->GetXaxis()->SetLabelSize(0.10);
  r->GetXaxis()->SetTitleOffset(1.2);
  r->Draw("e1");
  TLine *one = new TLine(r->GetXaxis()->GetXmin(), 1, r->GetXaxis()->GetXmax(), 1);
  one->SetLineStyle(2);
  one->Draw();

  c->SaveAs(outpng.c_str());
  delete c;
}

void drawMap2D(TH2D *hGen, TH2D *hProd, const std::string &prodLabel, const std::string &label,
               const std::string &outpng)
{
  if (!hGen || !hProd) return;
  TCanvas *c = new TCanvas("c_2d", "c_2d", 1400, 650);
  c->Divide(2, 1);
  gStyle->SetPalette(kBird);
  TH2D *hs[2] = {(TH2D *) hGen->Clone("gen2d"), (TH2D *) hProd->Clone("prod2d")};
  const char *titles[2] = {"genesis (luisvale)", prodLabel.c_str()};
  double zmax = 0;
  for (int i = 0; i < 2; i++)
  {
    // genesis's histogram was saved with its own TPaletteAxis, whose stored position would
    // override this canvas's margins -- drop it so both pads get a fresh, identical palette
    if (TObject *pal = hs[i]->GetListOfFunctions()->FindObject("palette")) hs[i]->GetListOfFunctions()->Remove(pal);
    if (hs[i]->Integral() > 0) hs[i]->Scale(1.0 / hs[i]->Integral());
    zmax = std::max(zmax, hs[i]->GetMaximum());
  }
  for (int i = 0; i < 2; i++)
  {
    c->cd(i + 1);
    gPad->SetTopMargin(0.10);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.25);
    gPad->SetLogz();
    hs[i]->SetTitle(";leading reco jet p_{T} [GeV];leading truth jet p_{T} [GeV]");
    hs[i]->SetMaximum(zmax);
    hs[i]->GetZaxis()->SetTitle("Normalized entries");
    hs[i]->GetZaxis()->SetTitleOffset(1.7);
    hs[i]->SetMinimum(1e-6);
    hs[i]->Draw("colz");
    TLine l;
    l.SetLineStyle(2);
    l.DrawLine(0, 0, 60, 60);
    TLatex lt;
    lt.SetNDC();
    lt.SetTextFont(42);
    lt.SetTextSize(0.030);
    TString l2(label);
    l2.ReplaceAll("\n", ", ");
    lt.DrawLatex(0.15, 0.93, (std::string(titles[i]) + ": " + l2.Data()).c_str());
  }
  c->SaveAs(outpng.c_str());
  delete c;
}

void compare_to_genesis(const std::string &treebase = "/sphenix/user/tmengel/cycleGAN/jewel_calo_trees",
                        const std::string &truthdir = "/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/qa_plots",
                        const std::string &outdir = "/sphenix/user/tmengel/cycleGAN/jewel_calo_trees/qa_plots/compare_genesis",
                        const std::string &genesisfile = "/sphenix/user/luisvale/projects/genesis/genesisObservables_R4.root",
                        const std::string &prodLabel = "this production (tmengel)")
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::AddDirectory(kFALSE);

  gSystem->mkdir(outdir.c_str(), kTRUE);

  TFile *fgen = TFile::Open(genesisfile.c_str(), "READ");
  if (!fgen || fgen->IsZombie())
  {
    std::cerr << "Cannot open genesis file " << genesisfile << std::endl;
    return;
  }

  for (const char *sim : {"vacuum", "norecoil"})
  {
    ProdHists p = buildProdHists(treebase, truthdir, sim);

    auto getGen = [&](const char *base) -> TH1 * { return (TH1 *) fgen->Get(Form("%s/%s_%s", sim, base, sim)); };

    const std::string simLabel = std::string("JEWEL ") + (std::string(sim) == "vacuum" ? "vacuum" : "no-recoil") + ", R = 0.4";
    const std::string fid = simLabel + "\nanti-k_{T}, |#eta^{jet}| < 0.7";
    const std::string win = fid + ", 20 < p_{T}^{jet} < 30 GeV";
    auto out = [&](const char *name) { return std::string(Form("%s/%s_%s.png", outdir.c_str(), sim, name)); };

    drawOverlay(getGen("hJetPtReco"), p.hJetPtReco, prodLabel, simLabel, out("jet_pt_reco"), true, true);
    drawOverlay(getGen("hJetPtTruth"), p.hJetPtTruth, prodLabel, simLabel, out("jet_pt_truth"), true, true);
    drawOverlay(getGen("hJetEta"), p.hJetEta, prodLabel, fid, out("jet_eta"), false);
    drawOverlay(getGen("hJetPhi"), p.hJetPhi, prodLabel, fid, out("jet_phi"), false);
    drawOverlay(getGen("hLeadJetPtReco"), p.hLeadJetPtReco, prodLabel, fid, out("lead_jet_pt_reco"), true, true);
    drawOverlay(getGen("hLeadJetPtTruth"), p.hLeadJetPtTruth, prodLabel, fid, out("lead_jet_pt_truth"), true, true);
    drawOverlay(getGen("hNJetsReco"), p.hNJetsReco, prodLabel, fid, out("n_jets_reco"), true);
    drawOverlay(getGen("hNJetsTruth"), p.hNJetsTruth, prodLabel, fid, out("n_jets_truth"), true);
    drawOverlay(getGen("hHepMCPartPt"), p.hHepMCPartPt, prodLabel, simLabel, out("gen_particle_pt"), true);
    drawOverlay(getGen("hHepMCPartEta"), p.hHepMCPartEta, prodLabel, simLabel, out("gen_particle_eta"), false);
    drawOverlay(getGen("hHepMCPartPhi"), p.hHepMCPartPhi, prodLabel, simLabel, out("gen_particle_phi"), false);
    drawMap2D((TH2D *) getGen("hLeadTruthVsRecoPt"), p.hLeadTruthVsRecoPt, prodLabel, fid, out("truth_vs_reco_leadpt"));

    if (p.hasConstituents)
    {
      drawOverlay(getGen("hGirthReco"), p.hGirthReco, prodLabel, win, out("girth_reco"), false);
      drawOverlay(getGen("hGirthTruth"), p.hGirthTruth, prodLabel, win, out("girth_truth"), false);
      drawOverlay(getGen("hNConstitReco"), p.hNConstitReco, prodLabel, win, out("n_constituents_reco"), false);
      drawOverlay(getGen("hNConstitTruth"), p.hNConstitTruth, prodLabel, win, out("n_constituents_truth"), false);
      drawOverlay(getGen("hGenPartPt"), p.hGenPartPt, prodLabel, fid, out("truthjet_constituent_pt"), true);
      drawOverlay(getGen("hGenPartEta"), p.hGenPartEta, prodLabel, fid, out("truthjet_constituent_eta"), false);
      drawOverlay(getGen("hGenPartPhi"), p.hGenPartPhi, prodLabel, fid, out("truthjet_constituent_phi"), false);

      // mean girth vs pT: profiles, drawn as markers without normalization
      TProfile *pg = (TProfile *) getGen("pGirthVsPt");
      if (pg)
      {
        TCanvas *c = new TCanvas("c_prof", "c_prof", 800, 650);
        gPad->SetLeftMargin(0.16);
        TProfile *a = (TProfile *) pg->Clone("gen_prof");
        TProfile *b = (TProfile *) p.pGirthVsPt->Clone("prod_prof");
        a->SetLineColor(kAzure + 2);
        a->SetMarkerColor(kAzure + 2);
        a->SetMarkerStyle(24);
        b->SetLineColor(kRed + 1);
        b->SetMarkerColor(kRed + 1);
        b->SetMarkerStyle(20);
        a->SetTitle(";reco jet p_{T} [GeV];#LT reco jet girth #GT");
        a->GetYaxis()->SetTitleOffset(1.5);
        a->SetMinimum(0);
        a->SetMaximum(0.30);
        a->GetXaxis()->SetRangeUser(10, 60);
        a->Draw("e1");
        b->Draw("e1 same");
        // data fall from top-left to bottom-right, so the bottom-left corner is empty
        TLegend *leg = new TLegend(0.20, 0.20, 0.56, 0.34);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(a, "genesis (luisvale)", "pe");
        leg->AddEntry(b, prodLabel.c_str(), "pe");
        leg->Draw();
        TLatex lt;
        lt.SetNDC();
        lt.SetTextFont(62);
        lt.SetTextSize(0.045);
        lt.DrawLatex(0.20, 0.86, "#it{#bf{sPHENIX}} Simulation");
        lt.SetTextFont(42);
        lt.SetTextSize(0.032);
        TObjArray *flines = TString(fid).Tokenize("\n");
        for (int k = 0; k < flines->GetEntries(); k++)
        {
          lt.DrawLatex(0.20, 0.80 - 0.05 * k, ((TObjString *) flines->At(k))->GetString().Data());
        }
        c->SaveAs(out("mean_girth_vs_pt").c_str());
        delete c;
      }
    }
  }

  std::cout << "Wrote comparison PNGs to " << outdir << std::endl;
}
