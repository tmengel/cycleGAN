// Jet-centred (Delta-eta, Delta-phi) energy-density maps for inclusive RECO jets, with the jet's
// towers split by calorimeter layer: EMCal, HCal (inner + outer), and their sum as a closure
// check. Two samples side by side, with their ratio and 1D projections.
//
// This is the calo-split companion to JetCenteredMaps.C. It is reco-only (truth particles have no
// calorimeter layer) and inclusive-only (no leading/subleading split), and it reproduces that
// macro's binning, normalisation, selection and panel layout, so the "all" component here is the
// same figure as JetCenteredMaps.C's jetmap_reco_inclusive.
//
// ---------------------------------------------------------------------------------------------
// Where the layer of a tower constituent comes from
//
// jettree's rjet_constit_src branch: the Jet::SRC each constituent was taken from, written out by
// JewelTreeWriter::FillJets straight from the jet's own comp_vec entry (comp.first). That is the
// container the jet reco actually read the tower out of, so it names the calorimeter with no
// reconstruction, no geometry and no assumptions -- which is why it is the only thing used here.
//
// The three sources the tower jets are clustered from are HCALIN_TOWERINFO (26),
// HCALOUT_TOWERINFO (27) and CEMC_TOWERINFO_RETOWER (28); see
// JewelTreeWriter::add_tower_constituent_source(). Any other value is counted and reported
// rather than silently folded into one of the components, so reconfiguring the tower sources
// shows up as a loud number instead of a wrong plot.
//
// Trees written before that branch existed cannot be split: the macro says so and stops. Do not
// try to reconstruct the layer from the kinematics -- tower eta is recomputed per event from the
// vertex (vz has an RMS of 64 cm here) and the three layers sit at different radii, so at a fixed
// eta bin the constituent eta spreads by ~0.4 and the layers overlap completely.

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
// ---------------------------------------------------------------------------------------------
// Normalisation: (1/N_jet) dpT / (dDelta-eta dDelta-phi) [GeV], divided by the jets entering the
// map and by the bin area. All three components share the same jets, hence the same denominator,
// so EMCal + HCal = all bin by bin. The projections integrate over the other axis, giving
// (1/N_jet) dpT / dDelta-eta [GeV].
//
// Jet selection (inclusive, as in JetCenteredMaps.C with inclusiveJets = true): every reco jet
// from jettree (already cut at pT > 10 GeV when written) with |eta| < kEtaFid and its OWN pT in
// [pt1min, pt1max). The event must have at least one fiducial truth jet, and with dijetVeto also
// a back-to-back truth dijet (|Delta-phi_12| >= 3pi/4) -- the same event selection the companion
// macro applies, so the two agree event for event.
//
// Output, per component (emcal, hcal, all): jetmap_reco_inclusive_<comp>.png (sample 1, sample 2,
// their ratio) and jetproj_reco_inclusive_<comp>.png (the Delta-eta and Delta-phi projections
// overlaid, each with a ratio panel). Then two figures about the split itself:
//   jetproj_calosplit_<sample>.png  one sample, the EMCal, HCal and total profiles overlaid, with
//                                   that sample's EMCal fraction underneath
//   jetproj_calofrac_ratio.png      both samples' EMCal and HCal fractions together, with the
//                                   sample ratio of each fraction underneath
//   jetfrac_average.png             the AVERAGE per-jet EMCal and HCal energy fractions: the
//                                   distribution of each jet's own share (left) and that share
//                                   against the jet's pT (right), both samples, each with the
//                                   sample ratio underneath. This is the per-jet average, a
//                                   different quantity from the per-bin density fractions the
//                                   two figures above show, and the means are also printed.
// and jet_centered_maps_calo.root with every histogram. The two fraction figures are drawn only
// over the bins holding at least kFracFloor of the peak tower density, since outside them there
// is no energy for a fraction to be a fraction of.
// Different selections overwrite each other, so scan pT bins by varying outdir.
//
// Usage:
//   root -l -b -q 'JetCenteredMapsCalo.C("<sample1 name>|<label>|<merged trees>",
//                                        "<sample2 ...>", "<outdir>", "<energy label>", nevents,
//                                        pt1min, pt1max, dijetVeto)'
// nevents <= 0 runs over all events. Sample 1 is the numerator of the ratio (vacuum). The sample
// spec takes a truthtree path as a fourth field for symmetry with JetCenteredMaps.C, but this
// macro does not read it.

#include <sPhenixStyle.C>

#include <TCanvas.h>
#include <TExec.h>
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
#include <TTree.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace
{
  // --- kept identical to JetCenteredMaps.C so the figures are directly comparable ---
  const double kWin = 0.6;      // half-width of the map (jets have R = 0.4)
  const double kEtaFid = 0.7;   // |eta| acceptance for the jet axis
  const int kNbins = 2.0 * (kWin / 0.09);   // calorimeter towers are ~0.09 wide
  const double kDphiMin = 0.75 * TMath::Pi();  // back-to-back requirement for the dijet veto

  const int kMapPalette = kInvertedDarkBodyRadiator;
  const int kRatioPalette = kTemperatureMap;   // diverging, so the ratio panel reads around 1
  const double kRatioMin = 0.0;                // symmetric about 1, which is the palette centre
  const double kRatioMax = 2.0;
  const double kRatioMaxRelErr = 0.5;          // 2D ratio bins noisier than this are left blank
  const double kHeadSize = 0.040;              // header text on the map panels

  // --- calo layers, as the Jet::SRC values rjet_constit_src stores ---
  enum Layer { kEMCal = 0, kIHCal = 1, kOHCal = 2, kNLayer = 3 };
  // Jet::SRC (jetbase/Jet.h) for each layer, in the same order
  const int kLayerSrc[kNLayer] = {28 /* CEMC_TOWERINFO_RETOWER */,
                                  26 /* HCALIN_TOWERINFO */,
                                  27 /* HCALOUT_TOWERINFO */};
  const int kUnknownSrc = -1;   // returned for a Jet::SRC that is none of the three

  const int kNcomp = 3;   // 0 = EMCal, 1 = HCal (inner + outer), 2 = all
  const char *kCompName[kNcomp] = {"emcal", "hcal", "all"};
  const char *kCompTitle[kNcomp] = {"EMCal towers", "HCal towers (in+out)", "all towers"};
  const char *kCompShort[kNcomp] = {"EMCal", "HCal", "all"};   // for the crowded fraction legend
  // the fraction figures only use bins holding at least this much of the peak tower density:
  // below it there is no energy for a fraction to be a fraction of
  const double kFracFloor = 1.0e-3;
  // which layers each component collects
  const bool kCompHas[kNcomp][kNLayer] = {{true, false, false},
                                          {false, true, true},
                                          {true, true, true}};

  struct Sample
  {
    std::string name, label, trees;
  };

  struct Selection
  {
    double pt1min, pt1max;
    bool dijetVeto;
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
    if (pt < sel.pt1min) return false;
    if (sel.pt1max > 0 && pt >= sel.pt1max) return false;
    return true;
  }

  // Event selection, on the truth jets, matching JetCenteredMaps.C's inclusive mode: at least one
  // fiducial truth jet, plus the dijet veto if asked. The pT window is a per-jet cut, applied to
  // the reco jets below, not an event veto.
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

  // Which calorimeter a tower constituent came from, straight off rjet_constit_src.
  int LayerOfSrc(int src)
  {
    for (int L = 0; L < kNLayer; L++)
    {
      if (src == kLayerSrc[L]) return L;
    }
    return kUnknownSrc;
  }

  // Integrates the density map over one axis: the projection sums densities, so multiplying by
  // the other axis' bin width turns that sum back into an integral.
  TH1D *Project(TH2D *h, bool alongX, const char *name, const std::string &normlab)
  {
    TH1D *p = alongX ? h->ProjectionX(name, 1, h->GetNbinsY())
                     : h->ProjectionY(name, 1, h->GetNbinsX());
    p->Scale(alongX ? h->GetYaxis()->GetBinWidth(1) : h->GetXaxis()->GetBinWidth(1));
    p->SetTitle(alongX
                    ? Form(";#Delta#eta;#LT (%s) dp_{T} / d#Delta#eta #GT [GeV]", normlab.c_str())
                    : Form(";#Delta#phi;#LT (%s) dp_{T} / d#Delta#phi #GT [GeV]", normlab.c_str()));
    return p;
  }

  std::string SelectionText(const Selection &sel)
  {
    std::string s = sel.pt1max > 0 ? Form("%0.1f < p_{T}^{jet,reco} < %.0f GeV", sel.pt1min, sel.pt1max)
                                   : Form("p_{T}^{jet,reco} > %.0f GeV", sel.pt1min);
    if (sel.dijetVeto) s += ", |#Delta#phi_{12}| > 3#pi/4";
    return s;
  }

  // Per-JET energy fractions, which is a different quantity from the per-bin density fractions
  // the map projections give: here each jet contributes one number, its own EMCal (or HCal) share
  // of its total tower pT, and those are then averaged over jets. That is the "average EMCal
  // energy fraction" of the jet sample, and it is what the medium modifies.
  struct FracHists
  {
    TH1D *dist[2] = {nullptr, nullptr};       // [0] EMCal share, [1] HCal share, per jet
    TProfile *vspt[2] = {nullptr, nullptr};   // the same two, averaged against the jet's own pT
    double sum[2] = {0, 0};                   // running mean, so the number can just be printed
    double sumsq[2] = {0, 0};
    long n = 0;
  };

  // per-layer tallies, reported per sample
  struct TagStats
  {
    long n[kNLayer] = {0, 0, 0};
    double pt[kNLayer] = {0, 0, 0};
    long nbad = 0;        // constituents whose Jet::SRC is none of the three layers
    double ptbad = 0;
  };
}  // namespace

// Fills the three component maps for one sample, and counts the jets that entered them.
// Returns false if the sample could not be read at all, so the caller stops instead of going on
// to draw and save a set of empty figures over the last good ones.
bool FillSampleCalo(const Sample &s, TH2D *maps[kNcomp], long &njet, long &nsel,
                    const Selection &sel, Long64_t nevents, TagStats &st, FracHists &fh)
{
  TFile ftrees(s.trees.c_str());
  TTree *jettree = (TTree *) ftrees.Get("jettree");
  if (!jettree)
  {
    printf("ERROR: %s has no jettree\n", s.trees.c_str());
    return false;
  }

  // The layer split is exactly this branch and nothing else. Trees written before it existed
  // cannot be split, and there is no reliable way to infer it from the kinematics, so stop.
  if (!jettree->GetBranch("rjet_constit_src"))
  {
    printf("ERROR: %s has no rjet_constit_src branch, so the tower jets cannot be split by\n"
           "       calorimeter. Regenerate the trees with a JewelTreeWriter that writes it.\n",
           s.trees.c_str());
    return false;
  }

  std::vector<float> *rjet_pt = nullptr, *rjet_eta = nullptr, *rjet_phi = nullptr;
  std::vector<float> *tjet_pt = nullptr, *tjet_eta = nullptr, *tjet_phi = nullptr;
  jettree->SetBranchAddress("rjet_pt", &rjet_pt);
  jettree->SetBranchAddress("rjet_eta", &rjet_eta);
  jettree->SetBranchAddress("rjet_phi", &rjet_phi);
  jettree->SetBranchAddress("tjet_pt", &tjet_pt);
  jettree->SetBranchAddress("tjet_eta", &tjet_eta);
  jettree->SetBranchAddress("tjet_phi", &tjet_phi);

  std::vector<std::vector<float>> *rc_pt = nullptr, *rc_eta = nullptr, *rc_phi = nullptr;
  std::vector<std::vector<int>> *rc_src = nullptr;
  jettree->SetBranchAddress("rjet_constit_pt", &rc_pt);
  jettree->SetBranchAddress("rjet_constit_eta", &rc_eta);
  jettree->SetBranchAddress("rjet_constit_phi", &rc_phi);
  jettree->SetBranchAddress("rjet_constit_src", &rc_src);

  Long64_t n = jettree->GetEntries();
  if (nevents > 0) n = std::min(n, nevents);

  nsel = 0;
  njet = 0;
  for (Long64_t i = 0; i < n; i++)
  {
    jettree->GetEntry(i);

    // the event selection is a truth-jet decision, the same one JetCenteredMaps.C takes
    const std::vector<size_t> tidx = RankFiducial(*tjet_pt, *tjet_eta);
    if (!PassEvent(*tjet_phi, tidx, sel)) continue;
    nsel++;

    const std::vector<size_t> ridx = RankFiducial(*rjet_pt, *rjet_eta);
    for (size_t q = 0; q < ridx.size(); q++)
    {
      const size_t k = ridx[q];
      if (!InPtWindow(rjet_pt->at(k), sel)) continue;
      njet++;

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

      // this jet's own pT split by layer, for its per-jet fractions below. Every constituent
      // counts, not just the ones inside the map window: the fraction is a property of the jet.
      double jptlay[kNLayer] = {0, 0, 0};

      for (size_t c = 0; c < cpt.size(); c++)
      {
        const int L = LayerOfSrc(csrc[c]);
        if (L == kUnknownSrc)
        {
          st.nbad++;
          st.ptbad += cpt[c];
          continue;
        }
        st.n[L]++;
        st.pt[L] += cpt[c];
        jptlay[L] += cpt[c];

        // the map window is applied after the layer lookup, so the tallies cover the whole jet
        const double dlong = ceta[c] - rjet_eta->at(k);
        if (std::fabs(dlong) > kWin) continue;
        const double dphi = TVector2::Phi_mpi_pi(cphi[c] - rjet_phi->at(k));
        if (std::fabs(dphi) > kWin) continue;

        for (int ic = 0; ic < kNcomp; ic++)
        {
          if (kCompHas[ic][L]) maps[ic]->Fill(dlong, dphi, cpt[c]);
        }
      }

      const double jptall = jptlay[kEMCal] + jptlay[kIHCal] + jptlay[kOHCal];
      if (jptall > 0)
      {
        const double f[2] = {jptlay[kEMCal] / jptall,
                             (jptlay[kIHCal] + jptlay[kOHCal]) / jptall};
        for (int c = 0; c < 2; c++)
        {
          fh.dist[c]->Fill(f[c]);
          fh.vspt[c]->Fill(rjet_pt->at(k), f[c]);
          fh.sum[c] += f[c];
          fh.sumsq[c] += f[c] * f[c];
        }
        fh.n++;
      }
    }
  }

  const double pttot = st.pt[kEMCal] + st.pt[kIHCal] + st.pt[kOHCal] + st.ptbad;
  printf("%s: %lld events read, %ld selected (%.1f%%); %ld inclusive reco jets in the pT window\n",
         s.name.c_str(), n, nsel, 100.0 * nsel / std::max(1LL, n), njet);
  printf("%s: tower constituents by layer -- EMCal %ld (%.2f%% of pT), IHCal %ld (%.2f%%), "
         "OHCal %ld (%.2f%%)\n",
         s.name.c_str(), st.n[kEMCal], 100.0 * st.pt[kEMCal] / std::max(1e-9, pttot),
         st.n[kIHCal], 100.0 * st.pt[kIHCal] / std::max(1e-9, pttot), st.n[kOHCal],
         100.0 * st.pt[kOHCal] / std::max(1e-9, pttot));
  if (st.nbad > 0)
  {
    printf("%s: WARNING %ld constituents (%.4f%% of pT) carry a Jet::SRC that is none of the "
           "three layers and were dropped -- the tower constituent sources this production was "
           "written with are not the ones kLayerSrc lists, so fix that before trusting the "
           "split\n",
           s.name.c_str(), st.nbad, 100.0 * st.ptbad / std::max(1e-9, pttot));
  }
  return true;
}

void JetCenteredMapsCalo(
  const std::string &spec1,
  const std::string &spec2,
  const std::string &outdir,
  const std::string &energyLabel = "JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV",
  Long64_t nevents = -1,
  double pt1min = 20,
  double pt1max = 30,
  bool dijetVeto = false
)
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kMapPalette);
  TH1::AddDirectory(kFALSE);
  gSystem->mkdir(outdir.c_str(), kTRUE);

  const Selection sel = {pt1min, pt1max, dijetVeto};
  const std::string seltext = SelectionText(sel);
  printf("event selection: %s\n", seltext.c_str());

  std::vector<Sample> samples = {ParseSample(spec1), ParseSample(spec2)};
  const std::string ratiolabel = samples[0].label + " / " + samples[1].label;
  const std::string shortratio = samples[0].name + " / " + samples[1].name;
  const std::string normlab = "1/N_{jet}";

  std::vector<std::array<TH2D *, kNcomp>> allmaps;
  std::vector<long> alln;
  std::vector<FracHists> allfrac;
  for (auto &s : samples)
  {
    TH2D *maps[kNcomp];
    for (int ic = 0; ic < kNcomp; ic++)
    {
      maps[ic] = new TH2D(Form("map_%s_reco_inclusive_%s", s.name.c_str(), kCompName[ic]),
                          ";#Delta#eta;#Delta#phi", kNbins, -kWin, kWin, kNbins, -kWin, kWin);
      maps[ic]->Sumw2();  // pT-weighted fills; the ratios need the errors
    }
    long njet = 0, nsel = 0;
    TagStats st;
    FracHists fh;
    // the jet pT axis spans the selection window, since that is all the jets there are
    const double ptlo = sel.pt1min;
    const double pthi = sel.pt1max > 0 ? sel.pt1max : sel.pt1min + 60.0;
    for (int c = 0; c < 2; c++)
    {
      fh.dist[c] = new TH1D(Form("fracdist_%s_%s", s.name.c_str(), kCompName[c]),
                            Form(";%s fraction of jet tower p_{T};jets (normalised)",
                                 kCompShort[c]),
                            50, 0.0, 1.0);
      fh.dist[c]->Sumw2();
      fh.vspt[c] = new TProfile(Form("fracvspt_%s_%s", s.name.c_str(), kCompName[c]),
                                Form(";p_{T}^{jet,reco} [GeV];#LT %s fraction #GT",
                                     kCompShort[c]),
                                10, ptlo, pthi);
    }
    if (!FillSampleCalo(s, maps, njet, nsel, sel, nevents, st, fh))
    {
      printf("nothing written\n");
      return;
    }
    std::array<TH2D *, kNcomp> m = {maps[0], maps[1], maps[2]};
    allmaps.push_back(m);
    alln.push_back(njet);
    allfrac.push_back(fh);

    // the headline number: the average share of the jet's tower pT each calorimeter carries
    for (int c = 0; c < 2; c++)
    {
      const double mean = fh.n > 0 ? fh.sum[c] / fh.n : 0.0;
      const double var = fh.n > 1 ? (fh.sumsq[c] / fh.n - mean * mean) * fh.n / (fh.n - 1) : 0.0;
      const double rms = var > 0 ? std::sqrt(var) : 0.0;
      printf("%s: average %-5s energy fraction per jet = %.4f +- %.4f (RMS %.4f, %ld jets)\n",
             s.name.c_str(), kCompShort[c], mean,
             fh.n > 0 ? rms / std::sqrt((double) fh.n) : 0.0, rms, fh.n);
    }
  }

  // Normalise up front, before anything is divided or projected. All three components come from
  // the same jets, so they share the denominator and stay additive: EMCal + HCal = all.
  for (size_t i = 0; i < samples.size(); i++)
  {
    for (int ic = 0; ic < kNcomp; ic++)
    {
      TH2D *h = allmaps[i][ic];
      if (alln[i] > 0)
      {
        h->Scale(1.0 / alln[i] / (h->GetXaxis()->GetBinWidth(1) * h->GetYaxis()->GetBinWidth(1)));
      }
    }
    printf("%s: every component normalised by its %ld jets\n", samples[i].name.c_str(), alln[i]);
  }

  TFile fout((outdir + "/jet_centered_maps_calo.root").c_str(), "RECREATE");

  // keep the Delta-eta projections of every component, per sample, for the overlay figure below
  std::vector<std::array<TH1D *, kNcomp>> profeta, profphi;

  for (int ic = 0; ic < kNcomp; ic++)
  {
    const std::string tag = std::string("reco_inclusive_") + kCompName[ic];

    // common z range across the two samples so they can be compared by eye
    double zmax = 0, zmin = 1e30;
    for (size_t i = 0; i < samples.size(); i++)
    {
      TH2D *h = allmaps[i][ic];
      zmax = std::max(zmax, h->GetMaximum());
      for (int bx = 1; bx <= h->GetNbinsX(); bx++)
      {
        for (int by = 1; by <= h->GetNbinsY(); by++)
        {
          const double v = h->GetBinContent(bx, by);
          if (v > 0) zmin = std::min(zmin, v);
        }
      }
    }
    zmin = std::max(zmin, zmax * 1e-5);

    TH2D *ratio = (TH2D *) allmaps[0][ic]->Clone(Form("ratio_%s", tag.c_str()));
    ratio->Divide(allmaps[1][ic]);  // an empty denominator gives 0, which COLZ leaves undrawn
    // A 2D panel shows no error bars, so the outer bins -- where both samples have a handful of
    // soft towers -- would otherwise read as real structure. Blank the ones that are noise.
    for (int bx = 1; bx <= ratio->GetNbinsX(); bx++)
    {
      for (int by = 1; by <= ratio->GetNbinsY(); by++)
      {
        const double v = ratio->GetBinContent(bx, by);
        if (v > 0 && ratio->GetBinError(bx, by) / v > kRatioMaxRelErr)
        {
          ratio->SetBinContent(bx, by, 0);
          ratio->SetBinError(bx, by, 0);
        }
      }
    }

    // ---- maps: sample 1, sample 2, ratio ----
    TCanvas c(Form("c_%s", tag.c_str()), "", 2100, 700);
    c.Divide(3, 1);
    for (int ipad = 0; ipad < 3; ipad++)
    {
      c.cd(ipad + 1);
      gPad->SetTopMargin(0.17);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.22);  // room for the palette labels AND the long z title
      gPad->SetBottomMargin(0.13);

      const bool isratio = ipad == 2;
      TH2D *h = isratio ? ratio : allmaps[ipad][ic];
      if (isratio)
      {
        h->SetMinimum(kRatioMin);
        h->SetMaximum(kRatioMax);
        h->GetZaxis()->SetTitle("Ratio");
      }
      else
      {
        h->SetMinimum(zmin);
        h->SetMaximum(zmax);
        h->GetZaxis()->SetTitle(
            Form("#LT (%s) dp_{T} / d#Delta#eta d#Delta#phi #GT [GeV]", normlab.c_str()));
      }
      h->GetZaxis()->SetTitleOffset(1.5);
      h->GetZaxis()->SetTitleSize(0.040);
      h->GetZaxis()->SetLabelSize(0.035);
      h->GetXaxis()->SetTitleOffset(1.1);
      h->GetYaxis()->SetTitleOffset(1.1);
      // The palette is global state in TStyle, so a pad can only claim its own by repainting:
      // the first draw builds the pad, the TExec swaps the palette, the second draw applies it.
      h->Draw("colz");
      (new TExec(Form("pal_%s_%d", tag.c_str(), ipad),
                 Form("gStyle->SetPalette(%d);", isratio ? kRatioPalette : kMapPalette)))
          ->Draw();
      h->Draw("colz same");

      TLatex lt;
      lt.SetNDC();
      lt.SetTextFont(42);
      lt.SetTextSize(kHeadSize);
      lt.DrawLatex(gPad->GetLeftMargin(), 0.955,
                   isratio ? ratiolabel.c_str() : samples[ipad].label.c_str());
      lt.DrawLatex(gPad->GetLeftMargin(), 0.910,
                   isratio ? Form("inclusive reco jet, %s", kCompTitle[ic])
                           : Form("inclusive reco jet, %s, %ld jets", kCompTitle[ic], alln[ipad]));
      lt.DrawLatex(gPad->GetLeftMargin(), 0.868, seltext.c_str());
    }
    c.SaveAs(Form("%s/jetmap_%s.png", outdir.c_str(), tag.c_str()));
    fout.cd();
    for (size_t i = 0; i < samples.size(); i++) allmaps[i][ic]->Write();
    ratio->Write();

    // ---- projections: Delta-eta and Delta-phi, both samples, each with a ratio panel ----
    TCanvas cp(Form("cp_%s", tag.c_str()), "", 1400, 800);
    const int col[2] = {kBlack, kRed + 1};
    const int mrk[2] = {20, 24};
    if (ic == 0)
    {
      profeta.resize(samples.size());
      profphi.resize(samples.size());
    }
    for (int v = 0; v < 2; v++)  // 0 = Delta-eta, 1 = Delta-phi
    {
      const char *vname = v == 0 ? "deta" : "dphi";
      const double x0 = 0.5 * v, x1 = 0.5 * (v + 1);
      cp.cd();
      TPad *pmain = new TPad(Form("main_%s_%s", tag.c_str(), vname), "", x0, 0.32, x1, 1.0);
      pmain->SetTopMargin(0.06);
      pmain->SetBottomMargin(0.02);
      pmain->SetLeftMargin(0.18);
      pmain->SetRightMargin(0.04);
      pmain->SetLogy();
      pmain->Draw();
      cp.cd();
      TPad *prat = new TPad(Form("rat_%s_%s", tag.c_str(), vname), "", x0, 0.0, x1, 0.32);
      prat->SetTopMargin(0.02);
      prat->SetBottomMargin(0.33);
      prat->SetLeftMargin(0.18);
      prat->SetRightMargin(0.04);
      prat->Draw();

      TH1D *p[2];
      double ymax = 0, ymin = 1e30;
      for (size_t i = 0; i < samples.size(); i++)
      {
        p[i] = Project(allmaps[i][ic], v == 0,
                       Form("proj_%s_%s_%s", samples[i].name.c_str(), tag.c_str(), vname), normlab);
        p[i]->SetLineColor(col[i]);
        p[i]->SetMarkerColor(col[i]);
        p[i]->SetMarkerStyle(mrk[i]);
        ymax = std::max(ymax, p[i]->GetMaximum());
        for (int b = 1; b <= p[i]->GetNbinsX(); b++)
        {
          if (p[i]->GetBinContent(b) > 0) ymin = std::min(ymin, p[i]->GetBinContent(b));
        }
        if (v == 0) profeta[i][ic] = p[i];
        else profphi[i][ic] = p[i];
      }
      ymin = std::max(ymin, ymax * 1e-5);

      pmain->cd();
      // the headroom keeps the peak clear of the header block, which log y makes easy to hit
      p[0]->GetYaxis()->SetRangeUser(0.5 * ymin, 100.0 * ymax);
      p[0]->GetXaxis()->SetLabelSize(0);
      p[0]->GetXaxis()->SetTitleSize(0);
      p[0]->GetYaxis()->SetTitleOffset(1.6);
      p[0]->Draw("e");
      p[1]->Draw("e same");

      TLegend *leg = new TLegend(0.64, 0.76, 0.96, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.040);
      for (size_t i = 0; i < samples.size(); i++) leg->AddEntry(p[i], samples[i].label.c_str(), "lp");
      leg->Draw();

      TLatex lt;
      lt.SetNDC();
      lt.SetTextFont(62);
      lt.SetTextSize(0.046);
      lt.DrawLatex(0.22, 0.920, "#it{#bf{sPHENIX}} Simulation");
      lt.SetTextFont(42);
      lt.SetTextSize(0.034);
      lt.DrawLatex(0.22, 0.875, Form("inclusive reco jet, %s", kCompTitle[ic]));
      lt.DrawLatex(0.22, 0.830, seltext.c_str());
      lt.DrawLatex(0.22, 0.785, energyLabel.c_str());

      prat->cd();
      TH1D *r = (TH1D *) p[0]->Clone(Form("projratio_%s_%s", tag.c_str(), vname));
      r->Divide(p[1]);
      r->GetYaxis()->SetRangeUser(kRatioMin, kRatioMax);
      r->GetYaxis()->SetTitle(shortratio.c_str());  // the full labels do not fit this pad
      r->GetYaxis()->SetNdivisions(505);
      // the ratio pad is 0.32 of the canvas against the main pad's 0.68, so its text grows to match
      const double sc = 0.68 / 0.32;
      r->GetXaxis()->SetLabelSize(p[0]->GetYaxis()->GetLabelSize() * sc);
      r->GetXaxis()->SetTitleSize(p[0]->GetYaxis()->GetTitleSize() * sc);
      r->GetXaxis()->SetTitleOffset(1.15);
      r->GetYaxis()->SetLabelSize(p[0]->GetYaxis()->GetLabelSize() * sc);
      r->GetYaxis()->SetTitleSize(p[0]->GetYaxis()->GetTitleSize() * sc * 0.8);
      r->GetYaxis()->SetTitleOffset(0.80);
      r->SetMarkerStyle(20);
      r->SetMarkerColor(kBlack);
      r->SetLineColor(kBlack);
      r->Draw("e");
      TLine *one = new TLine(-kWin, 1.0, kWin, 1.0);
      one->SetLineStyle(2);
      one->SetLineColor(kGray + 2);
      one->Draw();

      fout.cd();
      for (size_t i = 0; i < samples.size(); i++) p[i]->Write();
      r->Write();
    }
    cp.SaveAs(Form("%s/jetproj_%s.png", outdir.c_str(), tag.c_str()));
  }

  // ---- per sample: the EMCal, HCal and total profiles overlaid, with the EMCal fraction ----
  // This is the figure the split is for: where in the jet each calorimeter puts its energy.
  const int ccol[kNcomp] = {kAzure + 2, kOrange + 8, kBlack};
  const int cmrk[kNcomp] = {21, 22, 20};
  for (size_t i = 0; i < samples.size(); i++)
  {
    TCanvas cs(Form("cs_%s", samples[i].name.c_str()), "", 1400, 800);
    for (int v = 0; v < 2; v++)  // 0 = Delta-eta, 1 = Delta-phi
    {
      const char *vname = v == 0 ? "deta" : "dphi";
      const double x0 = 0.5 * v, x1 = 0.5 * (v + 1);
      cs.cd();
      TPad *pmain = new TPad(Form("csmain_%s_%s", samples[i].name.c_str(), vname), "", x0, 0.32, x1, 1.0);
      pmain->SetTopMargin(0.06);
      pmain->SetBottomMargin(0.02);
      pmain->SetLeftMargin(0.18);
      pmain->SetRightMargin(0.04);
      pmain->SetLogy();
      pmain->Draw();
      cs.cd();
      TPad *prat = new TPad(Form("csrat_%s_%s", samples[i].name.c_str(), vname), "", x0, 0.0, x1, 0.32);
      prat->SetTopMargin(0.02);
      prat->SetBottomMargin(0.33);
      prat->SetLeftMargin(0.18);
      prat->SetRightMargin(0.04);
      prat->Draw();

      TH1D *src[kNcomp];
      double ymax = 0, ymin = 1e30;
      for (int ic = 0; ic < kNcomp; ic++)
      {
        src[ic] = (TH1D *) (v == 0 ? profeta[i][ic] : profphi[i][ic])
                      ->Clone(Form("calosplit_%s_%s_%s", samples[i].name.c_str(), kCompName[ic], vname));
        src[ic]->SetLineColor(ccol[ic]);
        src[ic]->SetMarkerColor(ccol[ic]);
        src[ic]->SetMarkerStyle(cmrk[ic]);
        // the source projections were drawn once already, with an axis range stored on them, and
        // GetMaximum() hands back that stored value rather than the largest bin -- so clear it and
        // take the range off the bins
        src[ic]->SetMinimum(-1111);
        src[ic]->SetMaximum(-1111);
        for (int b = 1; b <= src[ic]->GetNbinsX(); b++)
        {
          const double y = src[ic]->GetBinContent(b);   // not v -- that is the axis selector
          if (y <= 0) continue;
          ymax = std::max(ymax, y);
          ymin = std::min(ymin, y);
        }
      }
      ymin = std::max(ymin, ymax * 1e-5);

      pmain->cd();
      // less headroom than the two-sample figures need: the three components span a much narrower
      // range than one sample's tails do, so 100x would leave four empty decades under the header
      src[kNcomp - 1]->GetYaxis()->SetRangeUser(0.5 * ymin, 20.0 * ymax);
      src[kNcomp - 1]->GetXaxis()->SetLabelSize(0);
      src[kNcomp - 1]->GetXaxis()->SetTitleSize(0);
      src[kNcomp - 1]->GetYaxis()->SetTitleOffset(1.6);
      src[kNcomp - 1]->Draw("e");
      for (int ic = 0; ic < kNcomp - 1; ic++) src[ic]->Draw("e same");

      // lower right, not the usual upper right: these profiles peak in the middle of the pad and
      // only span ~6 decades, so there is no headroom to put a legend over the peak without
      // either covering it or padding the axis out with empty decades
      TLegend *leg = new TLegend(0.60, 0.16, 0.95, 0.40);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.038);
      for (int ic = 0; ic < kNcomp; ic++) leg->AddEntry(src[ic], kCompTitle[ic], "lp");
      leg->Draw();

      TLatex lt;
      lt.SetNDC();
      lt.SetTextFont(62);
      lt.SetTextSize(0.046);
      lt.DrawLatex(0.22, 0.920, "#it{#bf{sPHENIX}} Simulation");
      lt.SetTextFont(42);
      lt.SetTextSize(0.034);
      lt.DrawLatex(0.22, 0.875, Form("%s, inclusive reco jet, %ld jets",
                                     samples[i].label.c_str(), alln[i]));
      lt.DrawLatex(0.22, 0.830, seltext.c_str());
      lt.DrawLatex(0.22, 0.785, energyLabel.c_str());

      // the EMCal fraction of the tower pT, bin by bin -- flat means the two calorimeters see the
      // same jet shape, a rise at large Delta-R means the shower spreads differently
      prat->cd();
      TH1D *frac = (TH1D *) src[kEMCal]->Clone(
          Form("emfrac_%s_%s", samples[i].name.c_str(), vname));
      frac->Divide(src[kNcomp - 1]);
      frac->GetYaxis()->SetRangeUser(0.0, 1.0);
      frac->GetYaxis()->SetTitle("EMCal / all");
      frac->GetYaxis()->SetNdivisions(505);
      const double sc = 0.68 / 0.32;
      frac->GetXaxis()->SetLabelSize(src[kNcomp - 1]->GetYaxis()->GetLabelSize() * sc);
      frac->GetXaxis()->SetTitleSize(src[kNcomp - 1]->GetYaxis()->GetTitleSize() * sc);
      frac->GetXaxis()->SetTitleOffset(1.15);
      frac->GetYaxis()->SetLabelSize(src[kNcomp - 1]->GetYaxis()->GetLabelSize() * sc);
      frac->GetYaxis()->SetTitleSize(src[kNcomp - 1]->GetYaxis()->GetTitleSize() * sc * 0.8);
      frac->GetYaxis()->SetTitleOffset(0.80);
      frac->SetMarkerStyle(20);
      frac->SetMarkerColor(kAzure + 2);
      frac->SetLineColor(kAzure + 2);
      frac->Draw("e");

      fout.cd();
      for (int ic = 0; ic < kNcomp; ic++) src[ic]->Write();
      frac->Write();
    }
    cs.SaveAs(Form("%s/jetproj_calosplit_%s.png", outdir.c_str(), samples[i].name.c_str()));
  }

  // ---- the EMCal and HCal fractions of the tower pT, both samples, with the sample ratio ----
  // Within one sample EMCal + HCal = all exactly, so the two fractions are mirror images about
  // 0.5 and carry the same information. Their ratios BETWEEN samples do not: EMCal_v/EMCal_n and
  // (1-EMCal_v)/(1-EMCal_n) respond differently to the same shift in the fraction, and the HCal
  // one is the more sensitive of the two wherever the EMCal dominates. Both are drawn.
  {
    const int nfrac = 2;   // 0 = EMCal / all, 1 = HCal / all; indexes the components directly
    const int fcol[nfrac] = {kAzure + 2, kOrange + 8};
    const int fmrk[nfrac][2] = {{21, 25}, {22, 26}};   // [component][sample]; filled = sample 1
    TCanvas cf("cf_calofrac", "", 1400, 800);
    for (int v = 0; v < 2; v++)  // 0 = Delta-eta, 1 = Delta-phi
    {
      const char *vname = v == 0 ? "deta" : "dphi";
      const char *vlab = v == 0 ? "#Delta#eta" : "#Delta#phi";
      const double x0 = 0.5 * v, x1 = 0.5 * (v + 1);
      cf.cd();
      TPad *pmain = new TPad(Form("cfmain_%s", vname), "", x0, 0.32, x1, 1.0);
      pmain->SetTopMargin(0.06);
      pmain->SetBottomMargin(0.02);
      pmain->SetLeftMargin(0.18);
      pmain->SetRightMargin(0.04);
      pmain->Draw();
      cf.cd();
      TPad *prat = new TPad(Form("cfrat_%s", vname), "", x0, 0.0, x1, 0.32);
      prat->SetTopMargin(0.02);
      prat->SetBottomMargin(0.33);
      prat->SetLeftMargin(0.18);
      prat->SetRightMargin(0.04);
      prat->Draw();

      // A fraction is meaningless where there is no energy to take a fraction of: the outermost
      // bins sit five orders below the peak density and their fractions scatter over the whole
      // [0, 1]. Keep the bins where either sample's total is within kFracFloor of its own peak,
      // and both draw and scale the fractions over that window only.
      TH1D *tot[2];
      for (size_t i = 0; i < samples.size(); i++) tot[i] = (v == 0 ? profeta[i][kNcomp - 1]
                                                                   : profphi[i][kNcomp - 1]);
      int blo = tot[0]->GetNbinsX() + 1, bhi = 0;
      for (size_t i = 0; i < samples.size(); i++)
      {
        const double peak = tot[i]->GetBinContent(tot[i]->GetMaximumBin());
        for (int b = 1; b <= tot[i]->GetNbinsX(); b++)
        {
          if (tot[i]->GetBinContent(b) < kFracFloor * peak) continue;
          blo = std::min(blo, b);
          bhi = std::max(bhi, b);
        }
      }

      // frac[c][i] = component c of sample i, divided by that sample's total. The numerator is a
      // subset of the denominator, so TH1::Divide's independent-errors assumption overestimates
      // these bars; that is the conservative direction and it is left alone.
      TH1D *frac[nfrac][2];
      double dev = 0;   // largest departure from 0.5 inside the window, for the y range
      for (int c = 0; c < nfrac; c++)
      {
        for (size_t i = 0; i < samples.size(); i++)
        {
          TH1D *num = (TH1D *) (v == 0 ? profeta[i][c] : profphi[i][c])
                          ->Clone(Form("calofrac_%s_%s_%s", samples[i].name.c_str(),
                                       kCompName[c], vname));
          num->SetMinimum(-1111);   // cleared for the same reason as in the overlay figure above
          num->SetMaximum(-1111);
          num->Divide(v == 0 ? profeta[i][kNcomp - 1] : profphi[i][kNcomp - 1]);
          num->SetLineColor(fcol[c]);
          num->SetMarkerColor(fcol[c]);
          num->SetMarkerStyle(fmrk[c][i]);
          num->SetTitle(Form(";%s;fraction of #LT dp_{T} / d%s #GT", vlab, vlab));
          frac[c][i] = num;
          for (int b = blo; b <= bhi; b++)
          {
            const double y = num->GetBinContent(b);
            if (y <= 0) continue;
            dev = std::max(dev, std::fabs(y - 0.5));
          }
        }
      }

      pmain->cd();
      // symmetric about 0.5, because the two fractions are mirror images about it, with enough
      // headroom for the legend to clear the EMCal peak
      const double half = std::min(0.5, std::max(0.12, 1.4 * dev));
      frac[0][0]->GetYaxis()->SetRangeUser(0.5 - half, 0.5 + half);
      frac[0][0]->GetXaxis()->SetRange(blo, bhi);
      frac[0][0]->GetXaxis()->SetLabelSize(0);
      frac[0][0]->GetXaxis()->SetTitleSize(0);
      frac[0][0]->GetYaxis()->SetTitleOffset(1.6);
      frac[0][0]->Draw("e");
      for (int c = 0; c < nfrac; c++)
      {
        for (size_t i = 0; i < samples.size(); i++)
        {
          if (c == 0 && i == 0) continue;
          frac[c][i]->Draw("e same");
        }
      }
      const double xlo = frac[0][0]->GetXaxis()->GetBinLowEdge(blo);
      const double xhi = frac[0][0]->GetXaxis()->GetBinUpEdge(bhi);
      TLine *half5 = new TLine(xlo, 0.5, xhi, 0.5);
      half5->SetLineStyle(2);
      half5->SetLineColor(kGray + 2);
      half5->Draw();

      // bottom right, not the usual top right: the header block reaches past the middle of the
      // pad, and the two fractions are mirror images about 0.5 so they fill the band the legend
      // would otherwise sit in. The short sample names are the ones the ratio axis already uses.
      TLegend *leg = new TLegend(0.58, 0.06, 0.96, 0.34);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.032);
      for (int c = 0; c < nfrac; c++)
      {
        for (size_t i = 0; i < samples.size(); i++)
        {
          leg->AddEntry(frac[c][i], Form("%s, %s", kCompShort[c], samples[i].name.c_str()), "lp");
        }
      }
      leg->Draw();

      TLatex lt;
      lt.SetNDC();
      lt.SetTextFont(62);
      lt.SetTextSize(0.046);
      lt.DrawLatex(0.22, 0.920, "#it{#bf{sPHENIX}} Simulation");
      lt.SetTextFont(42);
      lt.SetTextSize(0.034);
      lt.DrawLatex(0.22, 0.875, "inclusive reco jet, tower p_{T} fraction");
      lt.DrawLatex(0.22, 0.830, seltext.c_str());
      lt.DrawLatex(0.22, 0.785, energyLabel.c_str());

      // the ratio the fractions are here for: how much each calorimeter's share of the jet moves
      // between the two samples. Unlike the fractions themselves the two curves are independent.
      prat->cd();
      TH1D *fr[nfrac];
      double rdev = 0;
      for (int c = 0; c < nfrac; c++)
      {
        fr[c] = (TH1D *) frac[c][0]->Clone(Form("calofracratio_%s_%s", kCompName[c], vname));
        fr[c]->Divide(frac[c][1]);
        fr[c]->SetMarkerStyle(fmrk[c][0]);
        for (int b = blo; b <= bhi; b++)
        {
          const double y = fr[c]->GetBinContent(b);
          if (y <= 0) continue;
          rdev = std::max(rdev, std::fabs(y - 1.0));
        }
      }
      const double rhalf = std::min(1.0, std::max(0.04, 1.4 * rdev));
      fr[0]->GetYaxis()->SetRangeUser(1.0 - rhalf, 1.0 + rhalf);
      fr[0]->GetXaxis()->SetRange(blo, bhi);
      fr[0]->GetYaxis()->SetTitle(shortratio.c_str());  // the full labels do not fit this pad
      fr[0]->GetYaxis()->SetNdivisions(505);
      const double sc = 0.68 / 0.32;
      fr[0]->GetXaxis()->SetLabelSize(frac[0][0]->GetYaxis()->GetLabelSize() * sc);
      fr[0]->GetXaxis()->SetTitleSize(frac[0][0]->GetYaxis()->GetTitleSize() * sc);
      fr[0]->GetXaxis()->SetTitleOffset(1.15);
      fr[0]->GetYaxis()->SetLabelSize(frac[0][0]->GetYaxis()->GetLabelSize() * sc);
      fr[0]->GetYaxis()->SetTitleSize(frac[0][0]->GetYaxis()->GetTitleSize() * sc * 0.8);
      fr[0]->GetYaxis()->SetTitleOffset(0.80);
      fr[0]->Draw("e");
      fr[1]->Draw("e same");
      TLine *one = new TLine(-kWin, 1.0, kWin, 1.0);
      one->SetLineStyle(2);
      one->SetLineColor(kGray + 2);
      one->Draw();

      fout.cd();
      for (int c = 0; c < nfrac; c++)
      {
        for (size_t i = 0; i < samples.size(); i++) frac[c][i]->Write();
        fr[c]->Write();
      }
    }
    cf.SaveAs(Form("%s/jetproj_calofrac_ratio.png", outdir.c_str()));
  }

  // ---- the average per-jet EMCal and HCal energy fractions, and how the medium moves them ----
  // Left: the distribution of each jet's own fraction, both samples, means quoted. Right: the
  // same fractions averaged against the jet's pT. Both with the vacuum/norecoil ratio beneath.
  // EMCal + HCal = 1 jet by jet, so the two distributions are mirror images about 0.5; the point
  // of showing both is that their ratios are not mirror images.
  {
    const int fcol[2] = {kAzure + 2, kOrange + 8};
    const int fmrk[2][2] = {{21, 25}, {22, 26}};   // [component][sample]; filled = sample 1
    TCanvas ca("ca_avgfrac", "", 1400, 800);
    for (int v = 0; v < 2; v++)  // 0 = the distributions, 1 = versus jet pT
    {
      const char *vname = v == 0 ? "dist" : "vspt";
      const double x0 = 0.5 * v, x1 = 0.5 * (v + 1);
      ca.cd();
      TPad *pmain = new TPad(Form("camain_%s", vname), "", x0, 0.32, x1, 1.0);
      pmain->SetTopMargin(0.06);
      pmain->SetBottomMargin(0.02);
      pmain->SetLeftMargin(0.18);
      pmain->SetRightMargin(0.04);
      pmain->Draw();
      ca.cd();
      TPad *prat = new TPad(Form("carat_%s", vname), "", x0, 0.0, x1, 0.32);
      prat->SetTopMargin(0.02);
      prat->SetBottomMargin(0.33);
      prat->SetLeftMargin(0.18);
      prat->SetRightMargin(0.04);
      prat->Draw();

      TH1D *h[2][2];   // [component][sample]
      double ymax = 0;
      for (int c = 0; c < 2; c++)
      {
        for (size_t i = 0; i < samples.size(); i++)
        {
          if (v == 0)
          {
            // area-normalised, so two samples with different jet counts are comparable
            h[c][i] = (TH1D *) allfrac[i].dist[c]->Clone(
                Form("avgfrac_%s_%s_dist", samples[i].name.c_str(), kCompName[c]));
            if (h[c][i]->Integral() > 0) h[c][i]->Scale(1.0 / h[c][i]->Integral());
          }
          else
          {
            h[c][i] = allfrac[i].vspt[c]->ProjectionX(
                Form("avgfrac_%s_%s_vspt", samples[i].name.c_str(), kCompName[c]));
            h[c][i]->SetTitle(Form(";p_{T}^{jet,reco} [GeV];#LT energy fraction #GT"));
          }
          h[c][i]->SetLineColor(fcol[c]);
          h[c][i]->SetMarkerColor(fcol[c]);
          h[c][i]->SetMarkerStyle(fmrk[c][i]);
          for (int b = 1; b <= h[c][i]->GetNbinsX(); b++)
          {
            ymax = std::max(ymax, h[c][i]->GetBinContent(b) + h[c][i]->GetBinError(b));
          }
        }
      }

      pmain->cd();
      // the distributions start at zero so a linear scale from 0 is right; the pT profiles sit
      // around 0.4-0.6, so those get a window about 0.5 instead of a squashed band at mid-height
      if (v == 0) h[0][0]->GetYaxis()->SetRangeUser(0.0, 1.45 * ymax);
      else h[0][0]->GetYaxis()->SetRangeUser(0.30, 0.72);
      h[0][0]->GetXaxis()->SetLabelSize(0);
      h[0][0]->GetXaxis()->SetTitleSize(0);
      h[0][0]->GetYaxis()->SetTitleOffset(1.6);
      h[0][0]->Draw("e");
      for (int c = 0; c < 2; c++)
      {
        for (size_t i = 0; i < samples.size(); i++)
        {
          if (c == 0 && i == 0) continue;
          h[c][i]->Draw("e same");
        }
      }

      TLegend *leg = new TLegend(v == 0 ? 0.60 : 0.56, 0.62, 0.96, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.031);
      for (int c = 0; c < 2; c++)
      {
        for (size_t i = 0; i < samples.size(); i++)
        {
          const double m = allfrac[i].n > 0 ? allfrac[i].sum[c] / allfrac[i].n : 0.0;
          leg->AddEntry(h[c][i], Form("%s, %s: #LT f #GT = %.3f", kCompShort[c],
                                      samples[i].name.c_str(), m), "lp");
        }
      }
      leg->Draw();

      TLatex lt;
      lt.SetNDC();
      lt.SetTextFont(62);
      lt.SetTextSize(0.046);
      lt.DrawLatex(0.22, 0.920, "#it{#bf{sPHENIX}} Simulation");
      lt.SetTextFont(42);
      lt.SetTextSize(0.034);
      lt.DrawLatex(0.22, 0.875, "inclusive reco jet, per-jet energy fraction");
      lt.DrawLatex(0.22, 0.830, seltext.c_str());
      lt.DrawLatex(0.22, 0.785, energyLabel.c_str());

      prat->cd();
      TH1D *r[2];
      double rdev = 0;
      for (int c = 0; c < 2; c++)
      {
        r[c] = (TH1D *) h[c][0]->Clone(Form("avgfracratio_%s_%s", kCompName[c], vname));
        r[c]->Divide(h[c][1]);
        r[c]->SetMarkerStyle(fmrk[c][0]);
        for (int b = 1; b <= r[c]->GetNbinsX(); b++)
        {
          const double y = r[c]->GetBinContent(b);
          // ignore bins that are pure noise, which in the distribution tails is most of them
          if (y <= 0 || r[c]->GetBinError(b) > 0.10) continue;
          rdev = std::max(rdev, std::fabs(y - 1.0));
        }
      }
      const double rhalf = std::min(1.0, std::max(0.05, 1.4 * rdev));
      r[0]->GetYaxis()->SetRangeUser(1.0 - rhalf, 1.0 + rhalf);
      r[0]->GetYaxis()->SetTitle(shortratio.c_str());
      r[0]->GetYaxis()->SetNdivisions(505);
      const double sc = 0.68 / 0.32;
      r[0]->GetXaxis()->SetLabelSize(h[0][0]->GetYaxis()->GetLabelSize() * sc);
      r[0]->GetXaxis()->SetTitleSize(h[0][0]->GetYaxis()->GetTitleSize() * sc);
      r[0]->GetXaxis()->SetTitleOffset(1.15);
      r[0]->GetXaxis()->SetTitle(v == 0 ? "energy fraction of jet tower p_{T}"
                                        : "p_{T}^{jet,reco} [GeV]");
      r[0]->GetYaxis()->SetLabelSize(h[0][0]->GetYaxis()->GetLabelSize() * sc);
      r[0]->GetYaxis()->SetTitleSize(h[0][0]->GetYaxis()->GetTitleSize() * sc * 0.8);
      r[0]->GetYaxis()->SetTitleOffset(0.80);
      r[0]->Draw("e");
      r[1]->Draw("e same");
      TLine *one = new TLine(r[0]->GetXaxis()->GetXmin(), 1.0, r[0]->GetXaxis()->GetXmax(), 1.0);
      one->SetLineStyle(2);
      one->SetLineColor(kGray + 2);
      one->Draw();

      fout.cd();
      for (int c = 0; c < 2; c++)
      {
        for (size_t i = 0; i < samples.size(); i++) h[c][i]->Write();
        r[c]->Write();
      }
    }
    ca.SaveAs(Form("%s/jetfrac_average.png", outdir.c_str()));
  }

  // the modification, as a single number per calorimeter
  for (int c = 0; c < 2; c++)
  {
    const double m0 = allfrac[0].n > 0 ? allfrac[0].sum[c] / allfrac[0].n : 0.0;
    const double m1 = allfrac[1].n > 0 ? allfrac[1].sum[c] / allfrac[1].n : 0.0;
    // the two samples are independent, so the errors on the means add in quadrature
    double e[2] = {0, 0};
    for (int i = 0; i < 2; i++)
    {
      const FracHists &g = allfrac[i];
      const double mm = g.n > 0 ? g.sum[c] / g.n : 0.0;
      const double var = g.n > 1 ? (g.sumsq[c] / g.n - mm * mm) * g.n / (g.n - 1) : 0.0;
      e[i] = (g.n > 0 && var > 0) ? std::sqrt(var / g.n) : 0.0;
    }
    const double rat = m1 > 0 ? m0 / m1 : 0.0;
    const double erat = (m0 > 0 && m1 > 0)
                            ? rat * std::sqrt((e[0] / m0) * (e[0] / m0) + (e[1] / m1) * (e[1] / m1))
                            : 0.0;
    printf("average %-5s fraction: %s %.4f +- %.4f, %s %.4f +- %.4f, ratio %.4f +- %.4f "
           "(%+.2f%%)\n",
           kCompShort[c], samples[0].name.c_str(), m0, e[0], samples[1].name.c_str(), m1, e[1],
           rat, erat, 100.0 * (rat - 1.0));
  }

  // closure: the components were filled from the same jets with the same denominator, so the
  // EMCal and HCal maps have to add up to the total bin by bin
  for (size_t i = 0; i < samples.size(); i++)
  {
    double worst = 0;
    for (int bx = 1; bx <= allmaps[i][2]->GetNbinsX(); bx++)
    {
      for (int by = 1; by <= allmaps[i][2]->GetNbinsY(); by++)
      {
        const double tot = allmaps[i][2]->GetBinContent(bx, by);
        const double sum = allmaps[i][0]->GetBinContent(bx, by) +
                           allmaps[i][1]->GetBinContent(bx, by);
        if (tot > 0) worst = std::max(worst, std::fabs(sum - tot) / tot);
      }
    }
    printf("%s: EMCal + HCal vs all, worst relative bin difference %.2e\n",
           samples[i].name.c_str(), worst);
  }

  fout.Close();
  printf("wrote maps, projections and jet_centered_maps_calo.root to %s\n", outdir.c_str());
}
