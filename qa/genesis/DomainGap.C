// Quantifies how far apart the vacuum and no-recoil domains actually are, per jet pT bin, at
// both truth (generator particles, in rapidity) and reco (calorimeter towers) level.
//
// For feasibility the question is not "is there a difference in the mean" -- with 100k events
// every mean differs significantly -- but "how large is the difference compared to the jet-to-jet
// spread", because that is what a per-image translation has to learn. So each observable is
// reported as the two means, the relative shift, and Cohen's d = (m1-m2)/pooled sigma.
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TVector2.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
  const double kEtaFid = 0.7;
  const double kMatchDR = 0.3;   // reco-to-truth matching radius

  struct Acc
  {
    long n = 0;
    double s = 0, s2 = 0;
    void add(double x) { n++; s += x; s2 += x * x; }
    double mean() const { return n ? s / n : 0; }
    double sd() const { return n > 1 ? std::sqrt(std::max(0.0, (s2 - n * mean() * mean()) / (n - 1))) : 0; }
    double err() const { return n ? sd() / std::sqrt((double) n) : 0; }
  };

  // nconstit, girth (pT-weighted mean dR), core pT fraction (dR<0.1), total constituent pT
  struct Obs
  {
    Acc nc, girth, core, sumpt;
  };

  double Cohen(const Acc &a, const Acc &b)
  {
    if (a.n < 2 || b.n < 2) return 0;
    const double sp = std::sqrt(((a.n - 1) * a.sd() * a.sd() + (b.n - 1) * b.sd() * b.sd()) /
                                (a.n + b.n - 2));
    return sp > 0 ? (a.mean() - b.mean()) / sp : 0;
  }
}  // namespace

void Fill(const char *file, std::vector<Obs> &truth, std::vector<Obs> &reco,
          const std::vector<double> &edges, Long64_t nev, bool truthInRapidity)
{
  TFile f(file);
  TTree *t = (TTree *) f.Get("jettree");
  std::vector<float> *tpt = nullptr, *teta = nullptr, *tphi = nullptr;
  std::vector<float> *rpt = nullptr, *reta = nullptr, *rphi = nullptr;
  std::vector<std::vector<float>> *tcpt = nullptr, *tceta = nullptr, *tcphi = nullptr, *tce = nullptr;
  std::vector<std::vector<float>> *rcpt = nullptr, *rceta = nullptr, *rcphi = nullptr;
  t->SetBranchAddress("tjet_pt", &tpt); t->SetBranchAddress("tjet_eta", &teta); t->SetBranchAddress("tjet_phi", &tphi);
  t->SetBranchAddress("rjet_pt", &rpt); t->SetBranchAddress("rjet_eta", &reta); t->SetBranchAddress("rjet_phi", &rphi);
  t->SetBranchAddress("tjet_constit_pt", &tcpt); t->SetBranchAddress("tjet_constit_eta", &tceta);
  t->SetBranchAddress("tjet_constit_phi", &tcphi); t->SetBranchAddress("tjet_constit_e", &tce);
  t->SetBranchAddress("rjet_constit_pt", &rcpt); t->SetBranchAddress("rjet_constit_eta", &rceta);
  t->SetBranchAddress("rjet_constit_phi", &rcphi);

  Long64_t n = t->GetEntries();
  if (nev > 0) n = std::min(n, nev);

  long nsel = 0, nmatch = 0;
  for (Long64_t i = 0; i < n; i++)
  {
    t->GetEntry(i);

    // Reco-driven selection: take every fiducial reco jet in the pT window, pair it with the
    // nearest unclaimed truth jet within kMatchDR, and drop the pair if there is no partner.
    // Both levels then describe the same jets, so truth-vs-reco isolates the detector.
    std::vector<bool> used(tpt->size(), false);
    for (size_t r = 0; r < rpt->size(); r++)
    {
      if (std::fabs(reta->at(r)) >= kEtaFid) continue;
      int b = -1;
      for (size_t q = 0; q + 1 < edges.size(); q++)
        if (rpt->at(r) >= edges[q] && rpt->at(r) < edges[q + 1]) b = q;
      if (b < 0) continue;
      nsel++;

      int m = -1;
      double bestdr = kMatchDR;
      for (size_t u = 0; u < tpt->size(); u++)
      {
        if (used[u]) continue;
        const double de = teta->at(u) - reta->at(r);
        const double dp = TVector2::Phi_mpi_pi(tphi->at(u) - rphi->at(r));
        const double dr = std::sqrt(de * de + dp * dp);
        if (dr < bestdr) { bestdr = dr; m = (int) u; }
      }
      if (m < 0) continue;
      used[m] = true;
      nmatch++;

      for (int level = 0; level < 2; level++)
      {
        const size_t k = level == 0 ? (size_t) m : r;
        std::vector<float> *je = level == 0 ? teta : reta;
        std::vector<float> *jf = level == 0 ? tphi : rphi;
        auto *cp = level == 0 ? tcpt : rcpt;
        auto *ce = level == 0 ? tceta : rceta;
        auto *cf = level == 0 ? tcphi : rcphi;

        const std::vector<float> &a = cp->at(k);
        const std::vector<float> &e = ce->at(k);
        const std::vector<float> &p = cf->at(k);

        // The longitudinal coordinate. Reco towers are massless so eta IS y there; for truth it
        // is a choice: rapidity is the coordinate FastJet clustered in, pseudorapidity is the one
        // the detector measures. Comparing the two levels like for like means using eta for both.
        double axis = je->at(k);
        std::vector<double> lon(a.size());
        if (level == 0 && truthInRapidity)
        {
          TLorentzVector j4;
          std::vector<TLorentzVector> v(a.size());
          for (size_t c = 0; c < a.size(); c++)
          {
            v[c].SetPtEtaPhiE(a[c], e[c], p[c], tce->at(k)[c]);
            j4 += v[c];
          }
          axis = j4.Rapidity();
          for (size_t c = 0; c < a.size(); c++) lon[c] = v[c].Rapidity() - axis;
        }
        else
        {
          for (size_t c = 0; c < a.size(); c++) lon[c] = e[c] - axis;
        }

        double sum = 0, g = 0, core = 0;
        for (size_t c = 0; c < a.size(); c++)
        {
          const double dphi = TVector2::Phi_mpi_pi(p[c] - jf->at(k));
          const double dr = std::sqrt(lon[c] * lon[c] + dphi * dphi);
          sum += a[c];
          g += a[c] * dr;
          if (dr < 0.1) core += a[c];
        }
        if (sum <= 0) continue;
        Obs &o = level == 0 ? truth[b] : reco[b];
        o.nc.add((double) a.size());
        o.girth.add(g / sum);
        o.core.add(core / sum);
        o.sumpt.add(sum);
      }
    }
  }
  printf("%s: %ld reco jets selected, %ld matched (%.1f%%), %ld dropped\n",
         file, nsel, nmatch, 100.0 * nmatch / std::max(1L, nsel), nsel - nmatch);
}

void DomainGap(const char *fvac, const char *fnor, Long64_t nev = -1, bool truthInRapidity = false)
{
  const std::vector<double> edges = {10, 15, 20, 30, 40, 60};
  const size_t nb = edges.size() - 1;
  std::vector<Obs> tv(nb), rv(nb), tn(nb), rn(nb);
  Fill(fvac, tv, rv, edges, nev, truthInRapidity);
  Fill(fnor, tn, rn, edges, nev, truthInRapidity);
  printf("truth longitudinal coordinate: %s\n", truthInRapidity ? "rapidity (Delta-y)" : "pseudorapidity (Delta-eta)");

  const char *lvl[2] = {"TRUTH (generator particles)", "RECO (calorimeter towers)"};
  for (int level = 0; level < 2; level++)
  {
    printf("\n================ %s ================\n", lvl[level]);
    printf("%-11s %-18s %10s %10s %8s %8s\n", "pT bin", "observable", "vacuum", "no-recoil", "rel.diff", "Cohen d");
    for (size_t b = 0; b < nb; b++)
    {
      Obs &V = level == 0 ? tv[b] : rv[b];
      Obs &N = level == 0 ? tn[b] : rn[b];
      if (V.nc.n == 0 || N.nc.n == 0) continue;
      char bin[32];
      snprintf(bin, sizeof(bin), "%.0f-%.0f", edges[b], edges[b + 1]);
      struct Row { const char *name; Acc *a; Acc *c; };
      Row rows[4] = {{"N constituents", &V.nc, &N.nc},
                     {"girth <dR>", &V.girth, &N.girth},
                     {"core frac dR<0.1", &V.core, &N.core},
                     {"sum constit pT", &V.sumpt, &N.sumpt}};
      for (int r = 0; r < 4; r++)
      {
        const double m1 = rows[r].a->mean(), m2 = rows[r].c->mean();
        printf("%-11s %-18s %10.4f %10.4f %7.1f%% %8.3f\n", r == 0 ? bin : "", rows[r].name, m1, m2,
               100.0 * (m2 - m1) / m1, Cohen(*rows[r].a, *rows[r].c));
      }
      printf("%-11s %-18s %10ld %10ld\n", "", "jets", V.nc.n, N.nc.n);
    }
  }
}
