// Radial profiles of inclusive reco jets -- all towers, EMCal towers, HCal (inner + outer)
// towers -- for two samples, with their ratio.
//
// Every tower constituent of a selected jet is placed at its distance from the jet axis,
// DeltaR = sqrt(Delta-eta^2 + Delta-phi^2), and its pT histogrammed in DeltaR, in two
// normalisations:
//
//   <dpT/dDeltaR>                   [GeV]  the absolute pT flow, averaged over jets
//   rho(DeltaR) = <(1/pT) dpT/dDeltaR>     the jet shape: each jet normalised to its own pT
//
// each for all towers and split by layer. The layer is jettree's rjet_constit_src, the Jet::SRC
// the constituent was clustered from, so rho = rho_EMCal + rho_HCal and likewise for dpT/dDeltaR,
// bin by bin, by construction. The relative make-up of the jet shape is rho_EMCal / rho and
// rho_HCal / rho in each DeltaR bin; those two sum to 1.
//
// pT in rho is the scalar sum of the jet's constituent pT, not rjet_pt (the vector sum FastJet
// returns). With it the integral of rho over DeltaR is exactly 1 for every jet, so the three rho
// curves are fractions of the jet. The two differ by 0.7% on average (1.0066 +- 0.0035 over
// 3070 vacuum jets at 20-30 GeV).
//
// Errors are the standard error of the mean over JETS, from each jet's own per-bin values. Sumw2
// over constituents would treat the towers of one jet as independent and understate them. The
// EMCal fraction's error includes its covariance with the total it is part of; the HCal fraction
// is 1 minus it and has the same error. Sample ratios treat the two samples as independent.
//
// The binning is 0.05 in DeltaR over [0, 0.6]. That is half a tower, which is fine because the
// jet axis is continuous; nothing lands beyond 0.6 (max 0.572). The figures stop at the jet
// radius, kRjet = 0.4: the towers beyond it carry 0.05% of the jet pT -- they are clustered
// on tower centres, so a few sit just past the nominal radius -- and on a log axis those
// near-empty bins would stretch it by five decades and flatten everything else. They stay in
// the histograms, the ROOT file and the integral of rho, and the fraction beyond is printed.
//
// ---------------------------------------------------------------------------------------------
// Reading the vacuum / norecoil ratio: it is BELOW 1 at small DeltaR and ABOVE 1 at large
// DeltaR, i.e. the medium jets are NARROWER than the vacuum ones. That is the opposite of the
// naive quenching expectation and it is not a labelling or arithmetic error -- checked:
//   - the stored ratio really is vacuum divided by norecoil;
//   - norecoil really is the quenched sample (truth dijet imbalance <pT2/pT1> = 0.732 +- 0.003
//     against 0.819 +- 0.002 for vacuum);
//   - the narrowing is already there in the TRUTH jets, so it is not a detector effect.
//
// It is also not a pT-migration bias. "A medium jet seen at pT is really a vacuum jet from
// pT + delta" fails a two-observable test: the delta needed to explain the girth (+6 to +9 GeV)
// predicts MORE constituents, because vacuum jets gain ~0.14 particles per GeV, whereas norecoil
// jets have ~0.8 to ~1.0 FEWER particles than vacuum at the same pT.
//
// What the jets actually do, at fixed truth pT 25-30 GeV (both totalling 27.5 GeV):
//   norecoil has +2.1 GeV MORE pT inside DeltaR < 0.05 and less in every bin from 0.10 to 0.40,
//   and it is missing particles at every constituent pT, most of all the soft ones
//   (-0.43 per jet below 1 GeV, -0.25 from 1-2 GeV).
// So pT is not redistributed outward by the medium; the soft wide-angle part of the jet is
// missing, and what survives at a given pT is the hard core.
//
// That is what running JEWEL with the recoils discarded does: the energy the shower gives to the
// medium leaves with the scattering centres and never reappears in the cone, so the medium
// response that would broaden the jet is simply absent from the final state. The consequence for
// analysis: these samples cannot be used to study medium-induced broadening or soft wide-angle
// jet structure, and the sign of this ratio should not be quoted as a prediction for AuAu jet
// shapes. Confirming the mechanism directly would need a recoil-on sample to compare against.
//
// Selection as JetCenteredMaps(Calo).C in inclusive mode: events with at least one fiducial
// truth jet (and a back-to-back truth dijet, |Delta-phi_12| > 3pi/4, if dijetVeto); every reco jet
// with |eta| < 0.7 and its own pT in [ptmin, ptmax).
//
// Output in outdir:
//   jetshape_dpTdR.png   <dpT/dDeltaR> -- all, EMCal, HCal -- the two samples and their ratio
//   jetshape_rho.png     rho(DeltaR), same layout
//   jetshape_makeup.png  rho = rho_EMCal + rho_HCal stacked for each sample, and the fractions
//                        rho_EMCal/rho and rho_HCal/rho for both samples with their ratio
//   jetshape_calo.root   every histogram
//
// Usage:
//   root -l -b -q 'JetShapeCalo.C+("<name>|<label>|<merged trees>", "<sample 2 ...>",
//                                  "<outdir>", "<energy label>", nevents, ptmin, ptmax, dijetVeto)'
// nevents <= 0 runs over all events. Sample 1 is the ratio numerator. The trees must carry
// rjet_constit_src (trees_src/); against older trees the macro says so and writes nothing.

#include <sPhenixStyle.C>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace
{
  const double kEtaFid = 0.7;                  // |eta| acceptance for the jet axis
  const double kDphiMin = 0.75 * TMath::Pi();  // back-to-back requirement for the dijet veto
  const double kRmax = 0.6;
  const int kNr = 12;
  const double kRjet = 0.4;   // the anti-kT radius; the figures stop here

  enum Layer { kEMCal = 0, kIHCal = 1, kOHCal = 2, kNLayer = 3 };
  // Jet::SRC (jetbase/Jet.h) of each layer, in Layer order
  const int kLayerSrc[kNLayer] = {28 /* CEMC_TOWERINFO_RETOWER */, 26 /* HCALIN_TOWERINFO */,
                                  27 /* HCALOUT_TOWERINFO */};

  enum Comp { kAll = 0, kEm = 1, kHcal = 2, kNcomp = 3 };
  const char *kCompName[kNcomp] = {"all", "emcal", "hcal"};
  const char *kCompTitle[kNcomp] = {"all towers", "EMCal towers", "HCal towers (in+out)"};
  const bool kCompHas[kNcomp][kNLayer] = {{true, true, true},
                                          {true, false, false},
                                          {false, true, true}};

  enum Kind { kAbs = 0, kRho = 1, kNKind = 2 };

  struct Sample
  {
    std::string name, label, trees;
  };

  struct Selection
  {
    double ptmin, ptmax;
    bool dijetVeto;
  };

  // Per-jet sums over the selected jets, per normalisation, component and DeltaR bin: the jet's
  // value x, and x^2 for the spread. cross is x_EMCal * x_all, for the fraction's covariance.
  struct Acc
  {
    long njet = 0;
    double sum[kNKind][kNcomp][kNr] = {};
    double sumsq[kNKind][kNcomp][kNr] = {};
    double cross[kNKind][kNr] = {};
    double ptratio = 0;   // sum over jets of (scalar constituent pT sum) / rjet_pt
    double pttot = 0;     // tagged constituent pT, all jets
    double ptout = 0;     // ...of which beyond kRmax
    long nbad = 0;        // constituents whose Jet::SRC is none of the three layers
    double ptbad = 0;
  };

  Sample ParseSample(const std::string &spec)
  {
    TObjArray *f = TString(spec).Tokenize("|");
    Sample s;
    s.name = ((TObjString *) f->At(0))->GetString().Data();
    s.label = ((TObjString *) f->At(1))->GetString().Data();
    s.trees = ((TObjString *) f->At(2))->GetString().Data();
    return s;
  }

  // indices of the fiducial jets, sorted by descending pT
  std::vector<size_t> RankFiducial(const std::vector<float> &pt, const std::vector<float> &eta)
  {
    std::vector<size_t> idx;
    for (size_t k = 0; k < pt.size(); k++)
    {
      if (std::fabs(eta[k]) < kEtaFid) idx.push_back(k);
    }
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return pt[a] > pt[b]; });
    return idx;
  }

  bool InPtWindow(double pt, const Selection &sel)
  {
    if (pt < sel.ptmin) return false;
    if (sel.ptmax > 0 && pt >= sel.ptmax) return false;
    return true;
  }

  bool PassEvent(const std::vector<float> &phi, const std::vector<size_t> &idx,
                 const Selection &sel)
  {
    if (idx.empty()) return false;
    if (sel.dijetVeto)
    {
      if (idx.size() < 2) return false;
      if (std::fabs(TVector2::Phi_mpi_pi(phi[idx[1]] - phi[idx[0]])) < kDphiMin) return false;
    }
    return true;
  }

  int LayerOfSrc(int src)
  {
    for (int L = 0; L < kNLayer; L++)
    {
      if (src == kLayerSrc[L]) return L;
    }
    return -1;
  }

  std::string SelectionText(const Selection &sel)
  {
    std::string s = sel.ptmax > 0 ? Form("%0.0f < p_{T}^{jet,reco} < %.0f GeV", sel.ptmin, sel.ptmax)
                                  : Form("p_{T}^{jet,reco} > %.0f GeV", sel.ptmin);
    if (sel.dijetVeto) s += ", |#Delta#phi_{12}| > 3#pi/4";
    return s;
  }

  // mean over jets per DeltaR bin, divided by the bin width, with the standard error of the mean
  TH1D *MeanHist(const Acc &a, int kind, int ic, const std::string &name)
  {
    TH1D *h = new TH1D(name.c_str(), "", kNr, 0.0, kRmax);
    const double w = kRmax / kNr;
    const double n = a.njet;
    for (int b = 0; b < kNr && n > 1; b++)
    {
      const double m = a.sum[kind][ic][b] / n;
      const double var = (a.sumsq[kind][ic][b] / n - m * m) * n / (n - 1);
      h->SetBinContent(b + 1, m / w);
      h->SetBinError(b + 1, std::sqrt(std::max(0.0, var) / n) / w);
    }
    return h;
  }

  // rho_EMCal / rho per DeltaR bin, or rho_HCal / rho = 1 - that. The error propagates the
  // covariance of numerator and denominator, since the EMCal pT is part of the total.
  TH1D *FracHist(const Acc &a, bool hcal, const std::string &name)
  {
    TH1D *h = new TH1D(name.c_str(), "", kNr, 0.0, kRmax);
    const double n = a.njet;
    for (int b = 0; b < kNr && n > 1; b++)
    {
      const double me = a.sum[kRho][kEm][b] / n;
      const double ma = a.sum[kRho][kAll][b] / n;
      if (ma <= 0 || me <= 0) continue;
      const double cee = (a.sumsq[kRho][kEm][b] / n - me * me) * n / (n - 1);
      const double caa = (a.sumsq[kRho][kAll][b] / n - ma * ma) * n / (n - 1);
      const double cea = (a.cross[kRho][b] / n - me * ma) * n / (n - 1);
      const double f = me / ma;
      const double var = f * f * (cee / (me * me) + caa / (ma * ma) - 2.0 * cea / (me * ma)) / n;
      h->SetBinContent(b + 1, hcal ? 1.0 - f : f);
      h->SetBinError(b + 1, std::sqrt(std::max(0.0, var)));
    }
    return h;
  }

  void Style(TH1D *h, int col, int mrk)
  {
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(mrk);
  }

  void DrawHeader(const std::vector<std::string> &lines, double y0, double dy, double scale)
  {
    TLatex lt;
    lt.SetNDC();
    lt.SetTextFont(62);
    lt.SetTextSize(0.046 * scale);
    lt.DrawLatex(0.22, y0, "#it{#bf{sPHENIX}} Simulation");
    lt.SetTextFont(42);
    lt.SetTextSize(0.034 * scale);
    for (size_t i = 0; i < lines.size(); i++) lt.DrawLatex(0.22, y0 - dy * (i + 1), lines[i].c_str());
  }

  // One column of a canvas: the curves in a main pad and their ratios in a pad beneath it.
  void DrawColumn(TCanvas &c, int col, int ncol, const std::string &id,
                  const std::vector<TH1D *> &hs, const std::vector<std::string> &legs,
                  const std::vector<TH1D *> &rs, const std::string &ytitle, bool logy,
                  double ylo, double yhi, const std::string &rtitle,
                  const std::vector<std::string> &head, const double legbox[4], double xhi)
  {
    const double x0 = double(col) / ncol, x1 = double(col + 1) / ncol;
    c.cd();
    TPad *pm = new TPad(("m_" + id).c_str(), "", x0, 0.32, x1, 1.0);
    pm->SetTopMargin(0.06);
    pm->SetBottomMargin(0.02);
    pm->SetLeftMargin(0.18);
    pm->SetRightMargin(0.04);
    if (logy) pm->SetLogy();
    pm->Draw();
    c.cd();
    TPad *pr = new TPad(("r_" + id).c_str(), "", x0, 0.0, x1, 0.32);
    pr->SetTopMargin(0.02);
    pr->SetBottomMargin(0.33);
    pr->SetLeftMargin(0.18);
    pr->SetRightMargin(0.04);
    pr->Draw();

    pm->cd();
    TH1D *f = hs[0];
    f->SetTitle(Form(";#DeltaR;%s", ytitle.c_str()));
    f->GetXaxis()->SetRangeUser(0.0, xhi);
    f->GetYaxis()->SetRangeUser(ylo, yhi);
    f->GetXaxis()->SetLabelSize(0);
    f->GetXaxis()->SetTitleSize(0);
    f->GetYaxis()->SetTitleOffset(1.6);
    f->Draw("e");
    for (size_t i = 1; i < hs.size(); i++) hs[i]->Draw("e same");

    TLegend *leg = new TLegend(legbox[0], legbox[1], legbox[2], legbox[3]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.032);
    for (size_t i = 0; i < hs.size(); i++) leg->AddEntry(hs[i], legs[i].c_str(), "lp");
    leg->Draw();
    DrawHeader(head, 0.920, 0.045, 1.0);

    pr->cd();
    // symmetric about 1, scaled to the bins that are measured well enough to mean anything
    double dev = 0;
    for (TH1D *r : rs)
    {
      for (int b = 1; b <= r->FindBin(xhi - 1e-6); b++)
      {
        const double y = r->GetBinContent(b);
        if (y > 0 && r->GetBinError(b) < 0.10) dev = std::max(dev, std::fabs(y - 1.0));
      }
    }
    const double half = std::min(1.0, std::max(0.05, 1.4 * dev));
    TH1D *r0 = rs[0];
    const double sc = 0.68 / 0.32;   // the ratio pad is 0.32 of the canvas, the main pad 0.68
    r0->SetTitle(Form(";#DeltaR;%s", rtitle.c_str()));
    r0->GetXaxis()->SetRangeUser(0.0, xhi);
    r0->GetYaxis()->SetRangeUser(1.0 - half, 1.0 + half);
    r0->GetYaxis()->SetNdivisions(505);
    r0->GetXaxis()->SetLabelSize(f->GetYaxis()->GetLabelSize() * sc);
    r0->GetXaxis()->SetTitleSize(f->GetYaxis()->GetTitleSize() * sc);
    r0->GetXaxis()->SetTitleOffset(1.15);
    r0->GetYaxis()->SetLabelSize(f->GetYaxis()->GetLabelSize() * sc);
    r0->GetYaxis()->SetTitleSize(f->GetYaxis()->GetTitleSize() * sc * 0.8);
    r0->GetYaxis()->SetTitleOffset(0.80);
    r0->Draw("e");
    for (size_t i = 1; i < rs.size(); i++) rs[i]->Draw("e same");
    TLine *one = new TLine(0.0, 1.0, xhi, 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 2);
    one->Draw();
  }

  // One full-height column: rho_EMCal and rho_HCal stacked, with the total as markers on top --
  // which they reproduce by construction, so the markers sit exactly on the top of the stack.
  void DrawStack(TCanvas &c, int col, int ncol, const std::string &id, TH1D *em, TH1D *hc,
                 TH1D *tot, double yhi, const std::vector<std::string> &head, double xhi)
  {
    const double x0 = double(col) / ncol, x1 = double(col + 1) / ncol;
    c.cd();
    TPad *p = new TPad(("s_" + id).c_str(), "", x0, 0.0, x1, 1.0);
    p->SetTopMargin(0.04);
    p->SetBottomMargin(0.105);
    p->SetLeftMargin(0.18);
    p->SetRightMargin(0.04);
    p->Draw();
    p->cd();

    TH1D *e = (TH1D *) em->Clone(("st_em_" + id).c_str());
    TH1D *h = (TH1D *) hc->Clone(("st_hc_" + id).c_str());
    e->SetFillColor(kAzure + 2);
    e->SetLineColor(kAzure + 2);
    h->SetFillColor(kOrange + 8);
    h->SetLineColor(kOrange + 8);
    THStack *st = new THStack(("stack_" + id).c_str(),
                              ";#DeltaR;#rho(#DeltaR) = #LT (1/p_{T}) dp_{T} / d#DeltaR #GT");
    st->Add(e);   // EMCal at the bottom, HCal on top of it
    st->Add(h);
    st->SetMinimum(0.0);
    st->SetMaximum(yhi);
    st->Draw("hist");
    st->GetXaxis()->SetRangeUser(0.0, xhi);
    st->GetYaxis()->SetTitleOffset(1.6);
    // the full-height pad is 1/0.68 taller than the main pads, so its axis text is scaled down
    // to match their size on the page
    st->GetXaxis()->SetLabelSize(tot->GetYaxis()->GetLabelSize() * 0.68);
    st->GetXaxis()->SetTitleSize(tot->GetYaxis()->GetTitleSize() * 0.68);
    st->GetYaxis()->SetLabelSize(tot->GetYaxis()->GetLabelSize() * 0.68);
    st->GetYaxis()->SetTitleSize(tot->GetYaxis()->GetTitleSize() * 0.68);
    st->GetXaxis()->SetTitleOffset(1.1);
    TH1D *t = (TH1D *) tot->Clone(("st_tot_" + id).c_str());
    Style(t, kBlack, 20);
    t->Draw("e same");
    gPad->Modified();

    TLegend *leg = new TLegend(0.60, 0.72, 0.96, 0.84);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.032 * 0.68);
    leg->AddEntry(t, "all towers", "lp");
    leg->AddEntry(e, "EMCal towers", "f");
    leg->AddEntry(h, "HCal towers (in+out)", "f");
    leg->Draw();
    DrawHeader(head, 0.935, 0.031, 0.68);
  }
}  // namespace

// Accumulates the per-jet radial sums for one sample. False if the trees cannot be used.
bool FillShape(const Sample &s, const Selection &sel, Long64_t nevents, Acc &a)
{
  TFile ftrees(s.trees.c_str());
  TTree *t = (TTree *) ftrees.Get("jettree");
  if (!t)
  {
    printf("ERROR: %s has no jettree\n", s.trees.c_str());
    return false;
  }
  if (!t->GetBranch("rjet_constit_src"))
  {
    printf("ERROR: %s has no rjet_constit_src branch, so the tower jets cannot be split by\n"
           "       calorimeter. Use trees written with it (trees_src/).\n",
           s.trees.c_str());
    return false;
  }

  std::vector<float> *rjet_pt = nullptr, *rjet_eta = nullptr, *rjet_phi = nullptr;
  std::vector<float> *tjet_pt = nullptr, *tjet_eta = nullptr, *tjet_phi = nullptr;
  std::vector<std::vector<float>> *rc_pt = nullptr, *rc_eta = nullptr, *rc_phi = nullptr;
  std::vector<std::vector<int>> *rc_src = nullptr;
  t->SetBranchAddress("rjet_pt", &rjet_pt);
  t->SetBranchAddress("rjet_eta", &rjet_eta);
  t->SetBranchAddress("rjet_phi", &rjet_phi);
  t->SetBranchAddress("tjet_pt", &tjet_pt);
  t->SetBranchAddress("tjet_eta", &tjet_eta);
  t->SetBranchAddress("tjet_phi", &tjet_phi);
  t->SetBranchAddress("rjet_constit_pt", &rc_pt);
  t->SetBranchAddress("rjet_constit_eta", &rc_eta);
  t->SetBranchAddress("rjet_constit_phi", &rc_phi);
  t->SetBranchAddress("rjet_constit_src", &rc_src);

  Long64_t n = t->GetEntries();
  if (nevents > 0) n = std::min(n, nevents);
  const double w = kRmax / kNr;
  long nsel = 0;

  for (Long64_t i = 0; i < n; i++)
  {
    t->GetEntry(i);
    const std::vector<size_t> tidx = RankFiducial(*tjet_pt, *tjet_eta);
    if (!PassEvent(*tjet_phi, tidx, sel)) continue;
    nsel++;

    for (size_t k : RankFiducial(*rjet_pt, *rjet_eta))
    {
      if (!InPtWindow(rjet_pt->at(k), sel)) continue;
      const std::vector<float> &cpt = rc_pt->at(k);
      const std::vector<float> &ceta = rc_eta->at(k);
      const std::vector<float> &cphi = rc_phi->at(k);
      const std::vector<int> &csrc = rc_src->at(k);
      if (csrc.size() != cpt.size())
      {
        printf("ERROR: %s event %lld jet %zu has %zu constituents but %zu src entries\n",
               s.trees.c_str(), i, k, cpt.size(), csrc.size());
        return false;
      }

      // the jet's pT for rho: the scalar sum of the constituents that have a layer
      double ptj = 0;
      std::vector<int> lay(cpt.size());
      for (size_t c = 0; c < cpt.size(); c++)
      {
        lay[c] = LayerOfSrc(csrc[c]);
        if (lay[c] < 0)
        {
          a.nbad++;
          a.ptbad += cpt[c];
          continue;
        }
        ptj += cpt[c];
      }
      if (ptj <= 0) continue;
      a.njet++;
      a.ptratio += ptj / rjet_pt->at(k);

      double v[kNcomp][kNr] = {};   // this jet's pT per component and DeltaR bin
      for (size_t c = 0; c < cpt.size(); c++)
      {
        if (lay[c] < 0) continue;
        const double de = ceta[c] - rjet_eta->at(k);
        const double dp = TVector2::Phi_mpi_pi(cphi[c] - rjet_phi->at(k));
        const double dr = std::sqrt(de * de + dp * dp);
        a.pttot += cpt[c];
        const int b = (int) (dr / w);
        if (b >= kNr)
        {
          a.ptout += cpt[c];
          continue;
        }
        for (int ic = 0; ic < kNcomp; ic++)
        {
          if (kCompHas[ic][lay[c]]) v[ic][b] += cpt[c];
        }
      }

      for (int kind = 0; kind < kNKind; kind++)
      {
        const double scale = kind == kAbs ? 1.0 : 1.0 / ptj;
        for (int b = 0; b < kNr; b++)
        {
          for (int ic = 0; ic < kNcomp; ic++)
          {
            const double x = v[ic][b] * scale;
            a.sum[kind][ic][b] += x;
            a.sumsq[kind][ic][b] += x * x;
          }
          a.cross[kind][b] += (v[kEm][b] * scale) * (v[kAll][b] * scale);
        }
      }
    }
  }

  printf("%s: %lld events read, %ld selected; %ld inclusive reco jets in the pT window\n",
         s.name.c_str(), n, nsel, a.njet);
  printf("%s: <scalar constituent pT / rjet_pt> = %.4f; pT beyond DeltaR %.1f: %.4f%%\n",
         s.name.c_str(), a.njet > 0 ? a.ptratio / a.njet : 0.0, kRmax,
         100.0 * a.ptout / std::max(1e-9, a.pttot));
  if (a.nbad > 0)
  {
    printf("%s: WARNING %ld constituents (%.4f%% of pT) have a Jet::SRC that is none of the three "
           "layers and were left out\n",
           s.name.c_str(), a.nbad, 100.0 * a.ptbad / std::max(1e-9, a.pttot + a.ptbad));
  }
  return true;
}

void JetShapeCalo(
  const std::string &spec1,
  const std::string &spec2,
  const std::string &outdir,
  const std::string &energyLabel = "JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV",
  Long64_t nevents = -1,
  double ptmin = 20,
  double ptmax = 30,
  bool dijetVeto = false
)
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::AddDirectory(kFALSE);
  gSystem->mkdir(outdir.c_str(), kTRUE);

  const Selection sel = {ptmin, ptmax, dijetVeto};
  const std::string seltext = SelectionText(sel);
  printf("selection: %s\n", seltext.c_str());

  std::vector<Sample> samples = {ParseSample(spec1), ParseSample(spec2)};
  std::vector<Acc> acc(2);
  for (size_t i = 0; i < 2; i++)
  {
    if (!FillShape(samples[i], sel, nevents, acc[i]))
    {
      printf("nothing written\n");
      return;
    }
  }
  const std::string shortratio = samples[0].name + " / " + samples[1].name;

  // h[kind][comp][sample], frac[comp (0 EMCal, 1 HCal)][sample]
  TH1D *h[kNKind][kNcomp][2];
  TH1D *frac[2][2];
  const char *kindName[kNKind] = {"dpTdR", "rho"};
  for (int i = 0; i < 2; i++)
  {
    for (int kind = 0; kind < kNKind; kind++)
    {
      for (int ic = 0; ic < kNcomp; ic++)
      {
        h[kind][ic][i] = MeanHist(acc[i], kind, ic,
                                  Form("%s_%s_%s", kindName[kind], kCompName[ic],
                                       samples[i].name.c_str()));
      }
    }
    frac[0][i] = FracHist(acc[i], false, Form("frac_emcal_%s", samples[i].name.c_str()));
    frac[1][i] = FracHist(acc[i], true, Form("frac_hcal_%s", samples[i].name.c_str()));
  }

  // closures: EMCal + HCal = all bin by bin, and rho integrates to (1 - the pT beyond kRmax)
  for (int i = 0; i < 2; i++)
  {
    double worst = 0, integral = 0;
    for (int b = 1; b <= kNr; b++)
    {
      const double tot = h[kRho][kAll][i]->GetBinContent(b);
      const double sum = h[kRho][kEm][i]->GetBinContent(b) + h[kRho][kHcal][i]->GetBinContent(b);
      if (tot > 0) worst = std::max(worst, std::fabs(sum - tot) / tot);
      integral += tot * h[kRho][kAll][i]->GetBinWidth(b);
    }
    printf("%s: integral of rho = %.6f; EMCal + HCal vs all, worst relative bin difference %.1e\n",
           samples[i].name.c_str(), integral, worst);
  }

  // draw out to the jet radius and no further (see the header)
  const int blast = h[kAbs][kAll][0]->FindBin(kRjet - 1e-6);
  const double xhi = h[kAbs][kAll][0]->GetXaxis()->GetBinUpEdge(blast);
  for (int i = 0; i < 2; i++)
  {
    double beyond = 0;
    for (int b = blast + 1; b <= kNr; b++)
    {
      beyond += h[kRho][kAll][i]->GetBinContent(b) * h[kRho][kAll][i]->GetBinWidth(b);
    }
    printf("%s: fraction of the jet pT beyond DeltaR = %.1f (not drawn): %.4f%%\n",
           samples[i].name.c_str(), kRjet, 100.0 * beyond);
  }

  const int scol[2] = {kBlack, kRed + 1};
  const int smrk[2] = {20, 24};
  // bottom left: a falling profile on a log axis leaves that corner empty, while the header
  // fills the top left and the long sample labels do not fit beside it
  const double legLow[4] = {0.22, 0.05, 0.70, 0.19};
  TFile fout((outdir + "/jetshape_calo.root").c_str(), "RECREATE");

  // ---- <dpT/dDeltaR> and rho: all / EMCal / HCal, the two samples and their ratio ----
  for (int kind = 0; kind < kNKind; kind++)
  {
    const std::string ytitle = kind == kAbs ? "#LT dp_{T} / d#DeltaR #GT [GeV]"
                                            : "#rho(#DeltaR) = #LT (1/p_{T}) dp_{T} / d#DeltaR #GT";
    TCanvas c(Form("c_%s", kindName[kind]), "", 2100, 800);
    // a common y range across the three components, so EMCal and HCal read against the total
    double ymax = 0, ymin = 1e30;
    for (int ic = 0; ic < kNcomp; ic++)
    {
      for (int i = 0; i < 2; i++)
      {
        for (int b = 1; b <= blast; b++)
        {
          const double y = h[kind][ic][i]->GetBinContent(b);
          if (y <= 0) continue;
          ymax = std::max(ymax, y);
          ymin = std::min(ymin, y);
        }
      }
    }
    for (int ic = 0; ic < kNcomp; ic++)
    {
      for (int i = 0; i < 2; i++) Style(h[kind][ic][i], scol[i], smrk[i]);
      TH1D *r = (TH1D *) h[kind][ic][0]->Clone(Form("ratio_%s_%s", kindName[kind], kCompName[ic]));
      r->Divide(h[kind][ic][1]);
      Style(r, kBlack, 20);
      std::vector<std::string> legs;
      for (int i = 0; i < 2; i++) legs.push_back(Form("%s (%ld jets)", samples[i].label.c_str(),
                                                      acc[i].njet));
      DrawColumn(c, ic, kNcomp, Form("%s_%s", kindName[kind], kCompName[ic]),
                 {h[kind][ic][0], h[kind][ic][1]}, legs, {r}, ytitle, true, 0.3 * ymin,
                 8.0 * ymax, shortratio,
                 {Form("inclusive reco jet, %s", kCompTitle[ic]), seltext, energyLabel}, legLow, xhi);
      fout.cd();
      for (int i = 0; i < 2; i++) h[kind][ic][i]->Write();
      r->Write();
    }
    c.SaveAs(Form("%s/jetshape_%s.png", outdir.c_str(), kindName[kind]));
  }

  // ---- the make-up of rho: stacked per sample, then the two fractions compared ----
  {
    TCanvas c("c_makeup", "", 2100, 800);
    double ymax = 0;
    for (int i = 0; i < 2; i++)
    {
      for (int b = 1; b <= blast; b++) ymax = std::max(ymax, h[kRho][kAll][i]->GetBinContent(b));
    }
    for (int i = 0; i < 2; i++)
    {
      DrawStack(c, i, 3, samples[i].name, h[kRho][kEm][i], h[kRho][kHcal][i], h[kRho][kAll][i],
                1.45 * ymax,
                {Form("%s, inclusive reco jet, %ld jets", samples[i].label.c_str(), acc[i].njet),
                 seltext, energyLabel},
                xhi);
    }

    const int fcol[2] = {kAzure + 2, kOrange + 8};
    const int fmrk[2][2] = {{21, 25}, {22, 26}};   // [component][sample]; filled = sample 1
    std::vector<TH1D *> fs, rs;
    std::vector<std::string> legs;
    for (int f = 0; f < 2; f++)
    {
      for (int i = 0; i < 2; i++)
      {
        Style(frac[f][i], fcol[f], fmrk[f][i]);
        fs.push_back(frac[f][i]);
        legs.push_back(Form("%s, %s", f == 0 ? "#rho_{EMCal} / #rho" : "#rho_{HCal} / #rho",
                            samples[i].name.c_str()));
      }
      TH1D *r = (TH1D *) frac[f][0]->Clone(Form("ratio_frac_%s", f == 0 ? "emcal" : "hcal"));
      r->Divide(frac[f][1]);
      Style(r, fcol[f], fmrk[f][0]);
      rs.push_back(r);
    }
    // the fractions sit between ~0.3 and ~0.7 and mirror each other about 0.5, so the legend
    // goes in the empty band beneath them, on the right where both are near 0.5
    const double legFrac[4] = {0.56, 0.05, 0.96, 0.29};
    DrawColumn(c, 2, 3, "makeup_frac", fs, legs, rs, "fraction of #rho(#DeltaR)", false, 0.12,
               0.88, shortratio,
               {"inclusive reco jet, #rho = #rho_{EMCal} + #rho_{HCal}", seltext}, legFrac, xhi);
    fout.cd();
    for (TH1D *x : fs) x->Write();
    for (TH1D *x : rs) x->Write();
    c.SaveAs(Form("%s/jetshape_makeup.png", outdir.c_str()));
  }

  // ---- the numbers behind the figures ----
  printf("\n  DeltaR      rho: %-8s %-8s ratio        | rho_EMCal/rho: %-8s %-8s ratio\n",
         samples[0].name.c_str(), samples[1].name.c_str(), samples[0].name.c_str(),
         samples[1].name.c_str());
  for (int b = 1; b <= blast; b++)
  {
    const double r0 = h[kRho][kAll][0]->GetBinContent(b), r1 = h[kRho][kAll][1]->GetBinContent(b);
    const double e0 = h[kRho][kAll][0]->GetBinError(b), e1 = h[kRho][kAll][1]->GetBinError(b);
    const double f0 = frac[0][0]->GetBinContent(b), f1 = frac[0][1]->GetBinContent(b);
    const double g0 = frac[0][0]->GetBinError(b), g1 = frac[0][1]->GetBinError(b);
    const double rr = r1 > 0 ? r0 / r1 : 0, fr = f1 > 0 ? f0 / f1 : 0;
    const double err = (r0 > 0 && r1 > 0) ? rr * std::hypot(e0 / r0, e1 / r1) : 0;
    const double efr = (f0 > 0 && f1 > 0) ? fr * std::hypot(g0 / f0, g1 / f1) : 0;
    printf("  %.2f-%.2f  %8.4f %8.4f %.3f+-%.3f   |                %.4f   %.4f   %.3f+-%.3f\n",
           h[kRho][kAll][0]->GetBinLowEdge(b), h[kRho][kAll][0]->GetXaxis()->GetBinUpEdge(b), r0,
           r1, rr, err, f0, f1, fr, efr);
  }

  fout.Close();
  printf("wrote jetshape_dpTdR.png, jetshape_rho.png, jetshape_makeup.png and jetshape_calo.root "
         "to %s\n",
         outdir.c_str());
}
