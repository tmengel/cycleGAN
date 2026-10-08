// Truth-jet constituent multiplicity for two samples overlaid, in jet pT bins.
//
// For every fiducial truth jet (|eta| < 0.7, from jettree, already cut at pT > 10 GeV when
// written) the number of generator particles clustered into it -- the length of its
// tjet_constit_pt vector -- is histogrammed, one figure per jet pT bin. Inclusive jets, not
// just the leading one, so the pT bins are populated independently of the event's hardest jet.
//
// Each distribution is normalised to unit sum, so the y axis is the probability for a jet in
// that pT bin to have N constituents and the two samples are directly comparable.
//
// Mean and standard deviation are accumulated from the raw values, not from the binned
// histogram, so they are exact and unaffected by the axis range. Both are printed per bin and
// drawn on the figure, and the trend against jet pT is collected into a summary figure.
//
// Caveat: the genesis samples were filtered at the generator level to contain at least one jet
// with pT > 20 GeV and |eta| < 0.7, so the bins below 20 GeV are a biased selection of jets
// (they are the softer jets of events that had a hard one) and are not an inclusive jet sample.
// The bins at and above 20 GeV are the ones to compare.
//
// Usage:
//   root -l -b -q 'TruthConstituentMultiplicity.C+("<sample1 name>|<label>|<merged trees>|<truthtree>",
//                                                  "<sample2 ...>", "<outdir>", "<energy label>",
//                                                  nevents, "<pT bin edges>")'
// nevents <= 0 runs over all events. The bin edges are a comma-separated list, low to high.

#include "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.572/rootmacros/sPhenixStyle.C"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
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
  const double kEtaFid = 0.7;   // |eta| acceptance for the jet axis, as in JetCenteredMaps.C
  const double kMatchDR = 0.3;  // reco-to-truth jet matching radius, for matchToReco
  const int kNmax = 40;         // x-axis upper edge; the observed maximum is 37
  const double kFillAlpha = 0.35;

  struct Sample
  {
    std::string name, label, trees, truth;
  };

  // running sums, so the mean and sigma never depend on the histogram binning or range
  struct Stats
  {
    long n = 0;
    double sum = 0, sum2 = 0;
    double Mean() const { return n > 0 ? sum / n : 0; }
    double Sigma() const
    {
      return n > 1 ? std::sqrt(std::max(0.0, (sum2 - n * Mean() * Mean()) / (n - 1))) : 0;
    }
    double MeanErr() const { return n > 0 ? Sigma() / std::sqrt((double) n) : 0; }
    double SigmaErr() const { return n > 1 ? Sigma() / std::sqrt(2.0 * (n - 1)) : 0; }
  };

  Sample ParseSample(const std::string &spec)
  {
    TObjArray *f = TString(spec).Tokenize("|");
    Sample s;
    s.name = ((TObjString *) f->At(0))->GetString().Data();
    s.label = ((TObjString *) f->At(1))->GetString().Data();
    s.trees = ((TObjString *) f->At(2))->GetString().Data();
    s.truth = ((TObjString *) f->At(3))->GetString().Data();
    return s;
  }

  std::vector<double> ParseEdges(const std::string &spec)
  {
    std::vector<double> e;
    TObjArray *f = TString(spec).Tokenize(",");
    for (int i = 0; i < f->GetEntries(); i++)
    {
      e.push_back(((TObjString *) f->At(i))->GetString().Atof());
    }
    std::sort(e.begin(), e.end());
    return e;
  }
}  // namespace

void FillSample(const Sample &s, std::vector<TH1D *> &h, std::vector<Stats> &st,
                const std::vector<double> &edges, Long64_t nevents, bool matchToReco)
{
  TFile ftrees(s.trees.c_str());
  TTree *jettree = (TTree *) ftrees.Get("jettree");

  std::vector<float> *tjet_pt = nullptr, *tjet_eta = nullptr, *tjet_phi = nullptr;
  std::vector<float> *rjet_pt = nullptr, *rjet_eta = nullptr, *rjet_phi = nullptr;
  std::vector<std::vector<float>> *tc_pt = nullptr;
  jettree->SetBranchAddress("tjet_pt", &tjet_pt);
  jettree->SetBranchAddress("tjet_eta", &tjet_eta);
  jettree->SetBranchAddress("tjet_constit_pt", &tc_pt);
  if (matchToReco)
  {
    jettree->SetBranchAddress("tjet_phi", &tjet_phi);
    jettree->SetBranchAddress("rjet_pt", &rjet_pt);
    jettree->SetBranchAddress("rjet_eta", &rjet_eta);
    jettree->SetBranchAddress("rjet_phi", &rjet_phi);
  }
  long nsel = 0, nmatch = 0;

  Long64_t n = jettree->GetEntries();
  if (nevents > 0) n = std::min(n, nevents);

  for (Long64_t i = 0; i < n; i++)
  {
    jettree->GetEntry(i);

    if (matchToReco)
    {
      // Reco-driven: bin on the RECO jet's pT, then histogram the matched truth jet's
      // constituents, so the particle-level distribution describes the jets a detector-level
      // selection actually picks. Unmatched reco jets are dropped.
      std::vector<bool> used(tjet_pt->size(), false);
      for (size_t r = 0; r < rjet_pt->size(); r++)
      {
        if (std::fabs(rjet_eta->at(r)) >= kEtaFid) continue;
        int b = -1;
        for (size_t q = 0; q + 1 < edges.size(); q++)
          if (rjet_pt->at(r) >= edges[q] && rjet_pt->at(r) < edges[q + 1]) b = q;
        if (b < 0) continue;
        nsel++;
        int m = -1;
        double bestdr = kMatchDR;
        for (size_t u = 0; u < tjet_pt->size(); u++)
        {
          if (used[u]) continue;
          const double de = tjet_eta->at(u) - rjet_eta->at(r);
          const double dp = TVector2::Phi_mpi_pi(tjet_phi->at(u) - rjet_phi->at(r));
          const double dr = std::sqrt(de * de + dp * dp);
          if (dr < bestdr) { bestdr = dr; m = (int) u; }
        }
        if (m < 0) continue;
        used[m] = true;
        nmatch++;
        const double nc = tc_pt->at((size_t) m).size();
        h[b]->Fill(nc);
        st[b].n++;
        st[b].sum += nc;
        st[b].sum2 += nc * nc;
      }
      continue;
    }

    for (size_t k = 0; k < tjet_pt->size(); k++)
    {
      if (std::fabs(tjet_eta->at(k)) >= kEtaFid) continue;
      const double pt = tjet_pt->at(k);
      for (size_t b = 0; b + 1 < edges.size(); b++)
      {
        if (pt < edges[b] || pt >= edges[b + 1]) continue;
        const double nc = tc_pt->at(k).size();
        h[b]->Fill(nc);
        st[b].n++;
        st[b].sum += nc;
        st[b].sum2 += nc * nc;
      }
    }
  }
  if (matchToReco)
  {
    printf("%s: %lld events read; %ld reco jets selected, %ld matched (%.1f%%), %ld dropped\n",
           s.name.c_str(), n, nsel, nmatch, 100.0 * nmatch / std::max(1L, nsel), nsel - nmatch);
  }
  else
  {
    printf("%s: %lld events read\n", s.name.c_str(), n);
  }
}

void TruthConstituentMultiplicity(const std::string &spec1, const std::string &spec2,
                                  const std::string &outdir,
                                  const std::string &energyLabel = "JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV",
                                  Long64_t nevents = -1,
                                  const std::string &ptbins = "10,15,20,30,40,60",
                                  bool matchToReco = false)
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::AddDirectory(kFALSE);
  gSystem->mkdir(outdir.c_str(), kTRUE);

  const std::vector<Sample> samples = {ParseSample(spec1), ParseSample(spec2)};
  const std::vector<double> edges = ParseEdges(ptbins);
  const size_t nbin = edges.size() - 1;
  const int col[2] = {kBlue + 1, kRed + 1};

  std::vector<std::vector<TH1D *>> h(samples.size());
  std::vector<std::vector<Stats>> st(samples.size());
  for (size_t i = 0; i < samples.size(); i++)
  {
    st[i].resize(nbin);
    for (size_t b = 0; b < nbin; b++)
    {
      h[i].push_back(new TH1D(Form("nconstit_%s_pt%.0f_%.0f", samples[i].name.c_str(), edges[b], edges[b + 1]),
                              ";N_{constituents};Normalized Counts", kNmax + 1, -0.5, kNmax + 0.5));
    }
    FillSample(samples[i], h[i], st[i], edges, nevents, matchToReco);
  }

  printf("\ntruth jet constituent multiplicity, |eta_jet| < %.1f%s\n", kEtaFid,
         matchToReco ? ", binned on matched RECO jet pT" : "");
  printf("%-14s %-24s %8s %10s %10s\n", "pT bin [GeV]", "sample", "jets", "<N>", "sigma");
  for (size_t b = 0; b < nbin; b++)
  {
    for (size_t i = 0; i < samples.size(); i++)
    {
      printf("%5.0f - %-6.0f %-24s %8ld %6.2f +- %-4.2f %5.2f +- %.2f\n", edges[b], edges[b + 1],
             samples[i].name.c_str(), st[i][b].n, st[i][b].Mean(), st[i][b].MeanErr(),
             st[i][b].Sigma(), st[i][b].SigmaErr());
    }
  }

  TFile fout((outdir + "/truth_constituent_multiplicity.root").c_str(), "RECREATE");

  for (size_t b = 0; b < nbin; b++)
  {
    double ymax = 0;
    for (size_t i = 0; i < samples.size(); i++)
    {
      if (h[i][b]->Integral() > 0) h[i][b]->Scale(1.0 / h[i][b]->Integral());
      h[i][b]->SetLineColor(col[i]);
      h[i][b]->SetLineWidth(2);
      h[i][b]->SetFillColorAlpha(col[i], kFillAlpha);
      ymax = std::max(ymax, h[i][b]->GetMaximum());
    }

    TCanvas c(Form("c_pt%.0f_%.0f", edges[b], edges[b + 1]), "", 800, 650);
    gPad->SetTopMargin(0.12);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    h[0][b]->GetYaxis()->SetRangeUser(0, 1.45 * ymax);
    h[0][b]->GetXaxis()->SetTitleOffset(1.1);
    h[0][b]->GetYaxis()->SetTitleOffset(1.5);
    h[0][b]->Draw("hist");
    h[1][b]->Draw("hist same");

    TLatex lt;
    lt.SetNDC();
    lt.SetTextFont(42);
    lt.SetTextSize(0.037);
    lt.DrawLatex(gPad->GetLeftMargin(), 0.945,
                 Form("%s  vs.  %s", samples[0].label.c_str(), samples[1].label.c_str()));
    lt.DrawLatex(gPad->GetLeftMargin(), 0.900,
                 Form("anti-k_{T} R = 0.4  |#eta_{jet}| < %.1f  %s #in [%.0f, %.0f] GeV",
                      kEtaFid, matchToReco ? "p_{T}^{jet,reco}" : "p_{T}^{jet}",
                      edges[b], edges[b + 1]));

    TLegend leg(0.52, 0.72, 0.94, 0.86);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.040);
    for (size_t i = 0; i < samples.size(); i++) leg.AddEntry(h[i][b], samples[i].label.c_str(), "f");
    leg.Draw();

    lt.SetTextSize(0.034);
    for (size_t i = 0; i < samples.size(); i++)
    {
      lt.SetTextColor(col[i]);
      lt.DrawLatex(0.52, 0.66 - 0.05 * i,
                   Form("#LT N #GT = %.2f,  #sigma = %.2f", st[i][b].Mean(), st[i][b].Sigma()));
    }
    lt.SetTextColor(kBlack);
    lt.SetTextFont(62);
    lt.SetTextSize(0.040);
    lt.DrawLatex(0.52, 0.53, "#it{#bf{sPHENIX}} Simulation");
    lt.SetTextFont(42);
    lt.SetTextSize(0.034);
    lt.DrawLatex(0.52, 0.48, energyLabel.c_str());

    c.SaveAs(Form("%s/nconstit_truth_pt%.0f_%.0f.png", outdir.c_str(), edges[b], edges[b + 1]));
    fout.cd();
    for (size_t i = 0; i < samples.size(); i++) h[i][b]->Write();
  }

  // ---- summary: mean and sigma against jet pT ----
  TCanvas cs("c_summary", "", 1300, 600);
  cs.Divide(2, 1);
  TGraphErrors *g[2][2];  // [sample][0 = mean, 1 = sigma]
  for (int what = 0; what < 2; what++)
  {
    double lo = 1e30, hi = 0;
    for (size_t i = 0; i < samples.size(); i++)
    {
      g[i][what] = new TGraphErrors();
      g[i][what]->SetName(Form("%s_vs_pt_%s", what == 0 ? "mean" : "sigma", samples[i].name.c_str()));
      for (size_t b = 0; b < nbin; b++)
      {
        if (st[i][b].n == 0) continue;
        const double y = what == 0 ? st[i][b].Mean() : st[i][b].Sigma();
        const double ey = what == 0 ? st[i][b].MeanErr() : st[i][b].SigmaErr();
        const int p = g[i][what]->GetN();
        g[i][what]->SetPoint(p, 0.5 * (edges[b] + edges[b + 1]), y);
        g[i][what]->SetPointError(p, 0.5 * (edges[b + 1] - edges[b]), ey);
        lo = std::min(lo, y);
        hi = std::max(hi, y);
      }
      g[i][what]->SetLineColor(col[i]);
      g[i][what]->SetMarkerColor(col[i]);
      g[i][what]->SetMarkerStyle(i == 0 ? 20 : 24);
      g[i][what]->SetMarkerSize(1.3);
      fout.cd();
      g[i][what]->Write();
    }

    cs.cd(what + 1);
    gPad->SetTopMargin(0.10);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    const double pad = 0.35 * (hi - lo);
    TH1D *frame = new TH1D(Form("frame%d", what), "", 1, edges.front(), edges.back());
    frame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    frame->GetYaxis()->SetTitle(what == 0 ? "#LT N_{constituents} #GT"
                                          : "#sigma(N_{constituents})");
    frame->GetYaxis()->SetRangeUser(lo - pad, hi + 1.2 * pad);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->Draw();
    for (size_t i = 0; i < samples.size(); i++) g[i][what]->Draw("p same");

    TLegend *leg = new TLegend(0.50, 0.24, 0.93, 0.38);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    for (size_t i = 0; i < samples.size(); i++) leg->AddEntry(g[i][what], samples[i].label.c_str(), "lp");
    leg->Draw();

    TLatex lt;
    lt.SetNDC();
    lt.SetTextFont(62);
    lt.SetTextSize(0.045);
    lt.DrawLatex(0.20, 0.87, "#it{#bf{sPHENIX}} Simulation");
    lt.SetTextFont(42);
    lt.SetTextSize(0.036);
    lt.DrawLatex(0.20, 0.82, energyLabel.c_str());
    lt.DrawLatex(0.20, 0.77, Form("truth jets, anti-k_{T} R = 0.4, |#eta_{jet}| < %.1f", kEtaFid));
  }
  cs.SaveAs(Form("%s/nconstit_truth_summary.png", outdir.c_str()));

  fout.Close();
  printf("\nwrote figures and truth_constituent_multiplicity.root to %s\n", outdir.c_str());
}
