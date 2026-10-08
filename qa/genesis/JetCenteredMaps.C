// Jet-centred (Delta-eta, Delta-phi) energy-density maps for the leading and subleading jet,
// at truth and reco level, for two samples side by side, with their ratio and 1D projections.
//
// For each selected jet the map is filled with that jet's own constituents, weighted by pT:
//   reco  : the jet's calorimeter towers (rjet_constit_*)
//   truth : the jet's generator particles (tjet_constit_*)
// so the two levels are the same quantity, before and after the detector. Constituents carry
// their own eta/phi, so the maps need no tower-grid geometry.
//
// The maps do NOT fall exactly to zero at R = 0.4 at truth level. FastJet clusters in rapidity,
// but these branches store pseudorapidity, and for a massive particle at low pT eta is pushed
// away from zero relative to y. Soft baryons and kaons inside the jet therefore land at
// |Delta-eta| up to ~1.4. It is 0.1% of the jet pT, and it is truth-only: towers are massless,
// so y = eta exactly and the reco maps stop at ~0.5.
//
// useTruthRapidity fixes that by building the truth maps in Delta-y. It needs tjet_constit_e,
// which only the regenerated trees carry (trees_pid/merged/, not the v2 tag's trees/merged/);
// against the old trees the macro says so and stops. The jet's own y cannot be read back either
// -- it needs the jet mass, which is not stored -- so the axis is rebuilt by summing the
// constituent four-vectors, which is exactly what the jet is. Reco maps are unchanged, and
// Delta-phi is unaffected in both cases. Measured over 20k events: max |Delta-eta| 2.9 versus
// max |Delta-y| 0.50, and the constituents beyond 0.45 drop from 0.43% to 0.002%.
//
// Why not "all towers in a window around the jet": the towertree's etabin index does not
// convert to eta accurately enough to line up with the jet axis. Under the obvious uniform
// mapping (eta = -1.1 + (bin+0.5)*2.2/24) only 58% of the leading jet's pT lands within
// dR < 0.4 of its axis, and the resulting map is smeared flat in Delta-eta while staying sharp
// in Delta-phi. The phi mapping is fine. This does not affect the per-event images
// (draw_event_v2.C indexes them by bin, not by eta), but it does mean an eta calibration of
// that grid is needed before making window maps from towers.
//
// Normalisation: (1/N_jet) dpT / (dDelta-eta dDelta-phi) [GeV], i.e. divided by the number of
// jets entering the map and by the bin area, so the two samples and both ranks are comparable.
// The projections integrate over the other axis, giving (1/N_jet) dpT / dDelta-eta [GeV].
//
// With normToLeadTruth, every map of a sample is instead divided by that sample's accepted
// leading truth jets -- one per selected event, the same denominator for all four maps. Each
// map then reads per selected event rather than per jet, so a sample that finds fewer
// subleading or reco jets shows up in the ratio as less pT, instead of that difference being
// divided away. The per-map jet counts stay on the panels either way.
//
// Jet selection: from jettree (jets are already cut at pT > 10 GeV when written), |eta| < 0.7,
// ranked by pT. The subleading maps need at least two jets passing that cut.
//
// inclusiveJets drops the leading/subleading distinction: every fiducial jet whose OWN pT is in
// [pt1min, pt1max) enters one map per level, and the pT window stops being an event veto and
// becomes a per-jet cut (so reco jets are cut on reco pT, not on the leading truth jet). Output
// is 2 figures instead of 4, named *_inclusive. The dijet veto, if set, still applies as an
// event cut. normToLeadTruth then divides by the truth jets in the window rather than by one
// jet per event. Truth and reco jets are selected independently, each on its own pT: nothing
// pairs a reco jet with a truth jet, so the two maps carry their own jet counts.
//
// Event selection, decided on the TRUTH jets and applied to every map in the event, so the reco
// and subleading maps describe the same events as the truth leading map (lead/sublead mode only
// -- see inclusiveJets above):
//   - leading truth jet pT in [pt1min, pt1max)   (pt1max <= 0 means no upper edge)
//   - optionally (dijetVeto), a back-to-back truth dijet: the event is dropped unless there are
//     two fiducial truth jets with |Delta-phi_12| >= 3pi/4.
// Both use the |eta| < 0.7 fiducial list, the same one the maps are built from.
//
// Output, per level and rank: jetmap_*.png (sample 1, sample 2, and their ratio) and
// jetproj_*.png (the Delta-eta and Delta-phi projections overlaid, each with a ratio panel).
// Different selections overwrite each other, so scan pT bins by varying outdir.
//
// The truthtree lives in a different file from the trees; event i of one is event i of the
// other (verified: the leading truth-jet pT and the maximum generator-particle pT correlate at
// 0.56 event by event, versus 0.01 for any shifted pairing).
//
// Usage:
//   root -l -b -q 'JetCenteredMaps.C+("<sample1 name>|<label>|<merged trees>|<truthtree>",
//                                     "<sample2 ...>", "<outdir>", "<energy label>", nevents,
//                                     pt1min, pt1max, dijetVeto, normToLeadTruth,
//                                     useTruthRapidity, inclusiveJets)'
// nevents <= 0 runs over all events. Sample 1 is the numerator of the ratio (vacuum).

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
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace
{
  const double kWin = 0.6;        // half-width of the map (jets have R = 0.4)
  const double kEtaFid = 0.7;     // |eta| acceptance for the jet axis
  // const int kNbinsTruth = 20;     // fine binning for generator particles
  // const int kNbinsReco = 0.8/0.09;      // ~0.05; calorimeter towers are ~0.09 wide
  const int kNbinsReco = 2.0*(kWin/0.09);
  const int kNbinsTruth = 2.0*(kWin/0.09);
  const double kDphiMin = 0.75 * TMath::Pi();  // back-to-back requirement for the dijet veto

  const int kMapPalette = kInvertedDarkBodyRadiator;
  // const int kMapPalette = kBird;
  const int kRatioPalette = kTemperatureMap;   // diverging, so the ratio panel reads around 1
  const double kRatioMin = 0.0;                // symmetric about 1, which is the palette centre
  const double kRatioMax = 2.0;
  const double kRatioMaxRelErr = 0.5;          // 2D ratio bins noisier than this are left blank
  const double kHeadSize = 0.040;              // header text on the map panels

  struct Sample
  {
    std::string name, label, trees, truth;
  };

  // Event selection, evaluated on the truth jets only.
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
    s.truth = ((TObjString *) f->At(3))->GetString().Data();
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

  // applyLeadPt is false in inclusive mode, where the pT window is a per-jet cut instead of an
  // event veto on the leading truth jet. The dijet veto stays an event cut either way.
  bool PassEvent(const std::vector<float> &pt, const std::vector<float> &phi,
                 const std::vector<size_t> &idx, const Selection &sel, bool applyLeadPt)
  {
    if (idx.empty()) return false;
    if (applyLeadPt && !InPtWindow(pt[idx[0]], sel)) return false;
    if (sel.dijetVeto)
    {
      if (idx.size() < 2) return false;
      if (std::fabs(TVector2::Phi_mpi_pi(phi[idx[1]] - phi[idx[0]])) < kDphiMin) return false;
    }
    return true;
  }

  // "#Delta#eta" normally; "#Deltay" for the truth maps when they are built in rapidity
  const char *LongLabel(bool useY, int level) { return (useY && level == 0) ? "#Deltay" : "#Delta#eta"; }

  // Integrates the density map over one axis: the projection sums densities, so multiplying by
  // the other axis' bin width turns that sum back into an integral.
  TH1D *Project(TH2D *h, bool alongX, const char *name, const std::string &normlab,
                const char *longlab)
  {
    TH1D *p = alongX ? h->ProjectionX(name, 1, h->GetNbinsY())
                     : h->ProjectionY(name, 1, h->GetNbinsX());
    p->Scale(alongX ? h->GetYaxis()->GetBinWidth(1) : h->GetXaxis()->GetBinWidth(1));
    p->SetTitle(alongX ? Form(";%s;#LT (%s) dp_{T} / d%s #GT [GeV]", longlab, normlab.c_str(), longlab)
                       : Form(";#Delta#phi;#LT (%s) dp_{T} / d#Delta#phi #GT [GeV]", normlab.c_str()));
    return p;
  }

  std::string SelectionText(const Selection &sel, bool inclusive)
  {
    const char *sym = inclusive ? "p_{T}^{jet}" : "p_{T,1}^{truth}";
    std::string s = sel.pt1max > 0 ?
                      Form("%0.1f < %s < %.0f GeV", sel.pt1min, sym, sel.pt1max)
                        : Form("%s > %.0f GeV", sym, sel.pt1min);
    if (sel.dijetVeto) s += ", |#Delta#phi_{12}| > 3#pi/4";
    return s;
  }
}  // namespace

// Fills the four maps (truth/reco x leading/subleading) for one sample.
// njets[] counts the jets that entered each map, for the normalisation.
void FillSample(const Sample &s, TH2D *maps[4], long njets[4], long &nsel, const Selection &sel,
                Long64_t nevents, bool useY, bool inclusive)
{
  TFile ftrees(s.trees.c_str());
  TTree *jettree = (TTree *) ftrees.Get("jettree");

  std::vector<float> *rjet_pt = nullptr, *rjet_eta = nullptr, *rjet_phi = nullptr;
  std::vector<float> *tjet_pt = nullptr, *tjet_eta = nullptr, *tjet_phi = nullptr;
  jettree->SetBranchAddress("rjet_pt", &rjet_pt);
  jettree->SetBranchAddress("rjet_eta", &rjet_eta);
  jettree->SetBranchAddress("rjet_phi", &rjet_phi);
  jettree->SetBranchAddress("tjet_pt", &tjet_pt);
  jettree->SetBranchAddress("tjet_eta", &tjet_eta);
  jettree->SetBranchAddress("tjet_phi", &tjet_phi);

  std::vector<std::vector<float>> *rc_pt = nullptr, *rc_eta = nullptr, *rc_phi = nullptr;
  std::vector<std::vector<float>> *tc_pt = nullptr, *tc_eta = nullptr, *tc_phi = nullptr;
  jettree->SetBranchAddress("rjet_constit_pt", &rc_pt);
  jettree->SetBranchAddress("rjet_constit_eta", &rc_eta);
  jettree->SetBranchAddress("rjet_constit_phi", &rc_phi);
  jettree->SetBranchAddress("tjet_constit_pt", &tc_pt);
  jettree->SetBranchAddress("tjet_constit_eta", &tc_eta);
  jettree->SetBranchAddress("tjet_constit_phi", &tc_phi);

  // Only the regenerated trees (trees_pid/) carry the constituent energy rapidity needs.
  std::vector<std::vector<float>> *tc_e = nullptr;
  if (useY)
  {
    if (!jettree->GetBranch("tjet_constit_e"))
    {
      printf("ERROR: %s has no tjet_constit_e branch, so truth rapidity cannot be computed.\n"
             "       Point at the regenerated trees (trees_pid/merged/) or set useTruthRapidity = false.\n",
             s.trees.c_str());
      return;
    }
    jettree->SetBranchAddress("tjet_constit_e", &tc_e);
  }

  Long64_t n = jettree->GetEntries();
  if (nevents > 0) n = std::min(n, nevents);

  nsel = 0;
  for (Long64_t i = 0; i < n; i++)
  {
    jettree->GetEntry(i);

    // the event selection is a truth-jet decision, taken once and applied to every map below
    const std::vector<size_t> tidx = RankFiducial(*tjet_pt, *tjet_eta);
    if (!PassEvent(*tjet_pt, *tjet_phi, tidx, sel, !inclusive)) continue;
    nsel++;

    const std::vector<size_t> ridx = RankFiducial(*rjet_pt, *rjet_eta);

    for (int level = 0; level < 2; level++)  // 0 = truth, 1 = reco
    {
      const std::vector<float> &eta = level == 0 ? *tjet_eta : *rjet_eta;
      const std::vector<float> &phi = level == 0 ? *tjet_phi : *rjet_phi;
      const std::vector<size_t> &idx = level == 0 ? tidx : ridx;

      // Which jets go into which map slot. Inclusive mode puts every fiducial jet whose OWN pT
      // is in the window into slot 0 and leaves slot 1 empty; otherwise slots 0 and 1 are the
      // leading and subleading jet, with the pT window already applied as an event veto above.
      const std::vector<float> &jpt = level == 0 ? *tjet_pt : *rjet_pt;
      std::vector<std::pair<size_t, int>> todo;
      if (inclusive)
      {
        for (size_t q = 0; q < idx.size(); q++)
        {
          if (InPtWindow(jpt[idx[q]], sel)) todo.push_back(std::make_pair(idx[q], level * 2));
        }
      }
      else
      {
        for (int rank = 0; rank < 2; rank++)
        {
          if (idx.size() < (size_t) rank + 1) continue;
          todo.push_back(std::make_pair(idx[rank], level * 2 + rank));
        }
      }

      for (size_t d = 0; d < todo.size(); d++)
      {
        const size_t k = todo[d].first;
        const int im = todo[d].second;
        njets[im]++;

        const std::vector<float> &cpt = level == 0 ? tc_pt->at(k) : rc_pt->at(k);
        const std::vector<float> &ceta = level == 0 ? tc_eta->at(k) : rc_eta->at(k);
        const std::vector<float> &cphi = level == 0 ? tc_phi->at(k) : rc_phi->at(k);

        // In rapidity mode the truth jet axis has to be rebuilt: y needs the jet's mass, which
        // is not stored, but the jet is exactly the sum of its constituents. Reco is untouched --
        // towers are massless, so y and eta are the same thing there.
        const bool asY = useY && level == 0;
        double axisY = 0;
        if (asY)
        {
          TLorentzVector jet4;
          for (size_t c = 0; c < cpt.size(); c++)
          {
            TLorentzVector q;
            q.SetPtEtaPhiE(cpt[c], ceta[c], cphi[c], tc_e->at(k)[c]);
            jet4 += q;
          }
          axisY = jet4.Rapidity();
        }

        for (size_t c = 0; c < cpt.size(); c++)
        {
          double dlong;
          if (asY)
          {
            TLorentzVector q;
            q.SetPtEtaPhiE(cpt[c], ceta[c], cphi[c], tc_e->at(k)[c]);
            dlong = q.Rapidity() - axisY;
          }
          else
          {
            dlong = ceta[c] - eta[k];
          }
          if (std::fabs(dlong) > kWin) continue;
          const double dphi = TVector2::Phi_mpi_pi(cphi[c] - phi[k]);
          if (std::fabs(dphi) > kWin) continue;
          maps[im]->Fill(dlong, dphi, cpt[c]);
        }
      }
    }
  }
  if (inclusive)
  {
    printf("%s: %lld events read, %ld selected (%.1f%%); jets in the pT window: truth %ld, reco %ld\n",
           s.name.c_str(), n, nsel, 100.0 * nsel / n, njets[0], njets[2]);
  }
  else
  {
    printf("%s: %lld events read, %ld selected (%.1f%%); jets per map: truth lead %ld, truth sublead %ld, reco lead %ld, reco sublead %ld\n",
           s.name.c_str(), n, nsel, 100.0 * nsel / n, njets[0], njets[1], njets[2], njets[3]);
  }
}

void JetCenteredMaps(
  const std::string &spec1, 
  const std::string &spec2, 
  const std::string &outdir,
  const std::string &energyLabel = "JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV",
  Long64_t nevents = -1, 
  double pt1min = 20, 
  double pt1max = 30,
  bool dijetVeto = false, 
  bool normToLeadTruth = false,
  bool useTruthRapidity = false,
  bool inclusiveJets = false
)
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kMapPalette);
  TH1::AddDirectory(kFALSE);
  gSystem->mkdir(outdir.c_str(), kTRUE);

  const Selection sel = {pt1min, pt1max, dijetVeto};
  const std::string seltext = SelectionText(sel, inclusiveJets);
  printf("event selection: %s\n", seltext.c_str());

  std::vector<Sample> samples = {ParseSample(spec1), ParseSample(spec2)};
  const char *lvlname[2] = {"truth", "reco"};
  const char *lvltitle[2] = {"truth jet, particles", "reco jet, towers"};
  const char *rankname[2] = {inclusiveJets ? "inclusive" : "leading", "subleading"};
  // inclusive mode fills only slot 0 per level, so there is no subleading figure
  const int nrank = inclusiveJets ? 1 : 2;
  const std::string ratiolabel = samples[0].label + " / " + samples[1].label;
  const std::string shortratio = samples[0].name + " / " + samples[1].name;

  std::vector<std::array<TH2D *, 4>> allmaps;
  std::vector<std::array<long, 4>> alln;
  for (auto &s : samples)
  {
    TH2D *maps[4];
    long njets[4] = {0, 0, 0, 0};
    long nsel = 0;
    for (int level = 0; level < 2; level++)
    {
      for (int rank = 0; rank < 2; rank++)
      {
        const int nb = level == 0 ? kNbinsTruth : kNbinsReco;
        maps[level * 2 + rank] = new TH2D(Form("map_%s_%s_%s", s.name.c_str(), lvlname[level], rankname[rank]),
                                          Form(";%s;#Delta#phi", LongLabel(useTruthRapidity, level)),
                                          nb, -kWin, kWin, nb, -kWin, kWin);
        maps[level * 2 + rank]->Sumw2();  // pT-weighted fills; the ratios need the errors
      }
    }
    FillSample(s, maps, njets, nsel, sel, nevents, useTruthRapidity, inclusiveJets);
    std::array<TH2D *, 4> m = {maps[0], maps[1], maps[2], maps[3]};
    std::array<long, 4> nn = {njets[0], njets[1], njets[2], njets[3]};
    allmaps.push_back(m);
    alln.push_back(nn);
  }

  // Normalise every map up front, before anything is divided or projected. Each sample uses its
  // own denominator: per map, the jets in that map, or -- with normToLeadTruth -- the accepted
  // leading truth jets, i.e. one per selected event. The second makes every map a per-event
  // quantity on a common denominator, so a sample finding fewer subleading or reco jets shows up
  // in the ratio as less pT instead of being divided away.
  const std::string normlab = normToLeadTruth ? "1/N_{lead}^{truth}" : "1/N_{jet}";
  for (size_t i = 0; i < samples.size(); i++)
  {
    for (int im = 0; im < 4; im++)
    {
      TH2D *h = allmaps[i][im];
      const long nnorm = normToLeadTruth ? alln[i][0] : alln[i][im];
      if (nnorm > 0)
      {
        h->Scale(1.0 / nnorm / (h->GetXaxis()->GetBinWidth(1) * h->GetYaxis()->GetBinWidth(1)));
      }
    }
    if (normToLeadTruth)
    {
      printf("%s: every map normalised by %ld accepted leading truth jets\n",
             samples[i].name.c_str(), alln[i][0]);
    }
    else
    {
      printf("%s: each map normalised by its own jet count\n", samples[i].name.c_str());
    }
  }

  TFile fout((outdir + "/jet_centered_maps.root").c_str(), "RECREATE");
  for (int level = 0; level < 2; level++)
  {
    for (int rank = 0; rank < nrank; rank++)
    {
      const int im = level * 2 + rank;
      const std::string tag = std::string(lvlname[level]) + "_" + rankname[rank];

      // common z range across the two samples so they can be compared by eye
      double zmax = 0, zmin = 1e30;
      for (size_t i = 0; i < samples.size(); i++)
      {
        TH2D *h = allmaps[i][im];
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

      TH2D *ratio = (TH2D *) allmaps[0][im]->Clone(Form("ratio_%s", tag.c_str()));
      ratio->Divide(allmaps[1][im]);  // an empty denominator gives 0, which COLZ leaves undrawn
      // A 2D panel shows no error bars, so the outer bins -- where both samples have a handful of
      // soft constituents -- would otherwise read as real structure. Blank the ones that are noise.
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
        // fits the three header lines just above the frame, with room for the descenders of the
        // p_{T,1} subscripts, which otherwise clip into the top of the plot
        gPad->SetTopMargin(0.17);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.22);  // room for the palette labels AND the long z title
        gPad->SetBottomMargin(0.13);

        const bool isratio = ipad == 2;
        TH2D *h = isratio ? ratio : allmaps[ipad][im];
        if (isratio)
        {
          h->SetMinimum(kRatioMin);
          h->SetMaximum(kRatioMax);
          h->GetZaxis()->SetTitle("Ratio");
        }
        else
        {
          // gPad->SetLogz();
          h->SetMinimum(zmin);
          h->SetMaximum(zmax);
          h->GetZaxis()->SetTitle(
              Form("#LT (%s) dp_{T} / d%s d#Delta#phi #GT [GeV]", normlab.c_str(),
                   LongLabel(useTruthRapidity, level)));
        }
        // the z title has to clear the palette tick labels, which sPhenixStyle draws large
        h->GetZaxis()->SetTitleOffset(1.5);
        h->GetZaxis()->SetTitleSize(0.040);
        h->GetZaxis()->SetLabelSize(0.035);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.1);
        // The palette is global state in TStyle, so a pad can only claim its own by repainting:
        // the first draw builds the pad, the TExec swaps the palette, the second draw applies it.
        // Setting the palette before the single draw does not work -- every pad ends up the same.
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
                     isratio ? Form("%s %s", rankname[rank], lvltitle[level])
                             : Form("%s %s, %ld jets", rankname[rank], lvltitle[level], alln[ipad][im]));
        lt.DrawLatex(gPad->GetLeftMargin(), 0.868, seltext.c_str());
        // if (ipad == 0)
        // {
        //   lt.SetTextFont(62);
        //   lt.DrawLatex(gPad->GetLeftMargin(), 0.841, "#it{#bf{sPHENIX}} Simulation");
        //   lt.SetTextFont(42);
        //   lt.DrawLatex(gPad->GetLeftMargin(), 0.804, energyLabel.c_str());
        // }
      }
      c.SaveAs(Form("%s/jetmap_%s.png", outdir.c_str(), tag.c_str()));
      fout.cd();
      for (size_t i = 0; i < samples.size(); i++) allmaps[i][im]->Write();
      ratio->Write();

      // ---- projections: Delta-eta and Delta-phi, both samples, each with a ratio panel ----
      TCanvas cp(Form("cp_%s", tag.c_str()), "", 1400, 800);
      const int col[2] = {kBlack, kRed + 1};
      const int mrk[2] = {20, 24};
      for (int v = 0; v < 2; v++)  // 0 = Delta-eta, 1 = Delta-phi
      {
        const char *vname = v == 0 ? (useTruthRapidity && level == 0 ? "dy" : "deta") : "dphi";
        const double x0 = 0.5 * v, x1 = 0.5 * (v + 1);
        cp.cd();
        TPad *pmain = new TPad(Form("main_%s", vname), "", x0, 0.32, x1, 1.0);
        pmain->SetTopMargin(0.06);
        pmain->SetBottomMargin(0.02);
        pmain->SetLeftMargin(0.18);
        pmain->SetRightMargin(0.04);
        pmain->SetLogy();
        pmain->Draw();
        cp.cd();
        TPad *prat = new TPad(Form("rat_%s", vname), "", x0, 0.0, x1, 0.32);
        prat->SetTopMargin(0.02);
        prat->SetBottomMargin(0.33);
        prat->SetLeftMargin(0.18);
        prat->SetRightMargin(0.04);
        prat->Draw();

        TH1D *p[2];
        double ymax = 0, ymin = 1e30;
        for (size_t i = 0; i < samples.size(); i++)
        {
          p[i] = Project(allmaps[i][im], v == 0,
                         Form("proj_%s_%s_%s", samples[i].name.c_str(), tag.c_str(), vname), normlab,
                         LongLabel(useTruthRapidity, level));
          p[i]->SetLineColor(col[i]);
          p[i]->SetMarkerColor(col[i]);
          p[i]->SetMarkerStyle(mrk[i]);
          ymax = std::max(ymax, p[i]->GetMaximum());
          for (int b = 1; b <= p[i]->GetNbinsX(); b++)
          {
            if (p[i]->GetBinContent(b) > 0) ymin = std::min(ymin, p[i]->GetBinContent(b));
          }
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
        lt.DrawLatex(0.22, 0.875, Form("%s %s", rankname[rank], lvltitle[level]));
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
  }
  fout.Close();
  printf("wrote maps, projections and jet_centered_maps.root to %s\n", outdir.c_str());
}
