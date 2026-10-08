// compare_sim.C -- compare a new pass2 production (trees) with the v2 genesis trees.
//
//   root -b -q 'compare_sim.C+("<outdir>", "<new vac>", "<new med>", "<v2 vac>", "<v2 med>")'
//   (each input: one ROOT file or a wildcard such as ".../trees/vacuum/trees_*.root")
//
// Events are selected like the genesis filter, at truth level: leading truth jet (R=0.4) with
// pT > 20 GeV and |eta| < 0.7. All histograms are shapes normalised to the selected events.
// Three curves per sample: the new production, v2 (all vertices), and v2 with |vz| < 10 cm.
// The new samples have the vertex fixed at 0 while v2 has the 65 cm wide pp zero-angle vertex,
// so "v2, |vz|<10 cm" isolates everything except the vertex; the ratio panel is new / that.
//
// Writes <outdir>/compare_sim.pdf, page*.png and a summary table on stdout.

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace
{
  struct Obs
  {
    std::string name, title;
    int n;
    double lo, hi;
    bool logy;
  };

  const std::vector<Obs> kObs = {
      {"vz", ";vertex z [cm]", 60, -150, 150, true},
      {"tlead_pt", ";leading truth jet p_{T} [GeV]", 30, 20, 80, true},
      {"tlead_eta", ";leading truth jet #eta", 28, -0.7, 0.7, false},
      {"tlead_girth", ";leading truth jet girth", 40, 0, 0.2, false},
      {"tlead_nconst", ";leading truth jet constituents", 40, 0, 40, false},
      {"match_eff", ";leading truth jet has a reco match (#DeltaR<0.3)", 2, 0, 2, false},
      {"resp", ";p_{T}^{reco} / p_{T}^{truth} (matched leading jet)", 40, 0, 1.6, false},
      {"deta", ";#eta^{reco} - #eta^{truth} (matched leading jet)", 40, -0.2, 0.2, true},
      {"rmatch_nconst", ";matched reco jet constituents (towers)", 40, 0, 40, false},
      {"rlead_pt", ";leading reco jet p_{T} [GeV]", 40, 0, 80, true},
      {"e_em", ";EMCal energy per event [GeV]", 50, 0, 100, false},
      {"e_ih", ";iHCal energy per event [GeV]", 50, 0, 20, false},
      {"e_oh", ";oHCal energy per event [GeV]", 50, 0, 50, false},
      {"e_tot", ";total calorimeter energy per event [GeV]", 50, 0, 150, false},
      {"ntow", ";towers with E > 50 MeV", 50, 0, 250, false},
      {"eprof", ";tower #eta bin (image row);mean energy per event [GeV]", 24, 0, 24, false},
  };

  struct Sample
  {
    std::map<std::string, TH1D *> h;
    long nread{0}, nsel{0};
  };

  Sample fillSample(const std::string &file, const std::string &tag, double vzcut)
  {
    Sample s;
    for (const auto &o : kObs)
    {
      s.h[o.name] = new TH1D((o.name + "_" + tag).c_str(), o.title.c_str(), o.n, o.lo, o.hi);
      s.h[o.name]->SetDirectory(nullptr);
      s.h[o.name]->Sumw2();
    }
    // file may be a single ROOT file or a wildcard over per-job tree files (TChain)
    TChain *jt = new TChain("jettree"), *tt = new TChain("towertree"), *gt = new TChain("globaltree");
    jt->Add(file.c_str());
    tt->Add(file.c_str());
    gt->Add(file.c_str());
    float vz = 0;
    std::vector<float> *tpt = nullptr, *teta = nullptr, *tphi = nullptr, *rpt = nullptr, *reta = nullptr, *rphi = nullptr;
    std::vector<std::vector<float>> *tcpt = nullptr, *tceta = nullptr, *tcphi = nullptr, *rcpt = nullptr;
    std::vector<float> *em = nullptr, *ih = nullptr, *oh = nullptr;
    std::vector<int> *etabin = nullptr;
    gt->SetBranchAddress("vz", &vz);
    jt->SetBranchAddress("tjet_pt", &tpt);
    jt->SetBranchAddress("tjet_eta", &teta);
    jt->SetBranchAddress("tjet_phi", &tphi);
    jt->SetBranchAddress("tjet_constit_pt", &tcpt);
    jt->SetBranchAddress("tjet_constit_eta", &tceta);
    jt->SetBranchAddress("tjet_constit_phi", &tcphi);
    jt->SetBranchAddress("rjet_pt", &rpt);
    jt->SetBranchAddress("rjet_eta", &reta);
    jt->SetBranchAddress("rjet_phi", &rphi);
    jt->SetBranchAddress("rjet_constit_pt", &rcpt);
    tt->SetBranchAddress("emcalenergy", &em);
    tt->SetBranchAddress("ihcalenergy", &ih);
    tt->SetBranchAddress("ohcalenergy", &oh);
    tt->SetBranchAddress("etabin", &etabin);

    const long n = jt->GetEntries();
    for (long i = 0; i < n; ++i)
    {
      gt->GetEntry(i);
      jt->GetEntry(i);
      s.nread++;
      if (vzcut > 0 && std::fabs(vz) >= vzcut) continue;
      int it = -1;
      for (size_t k = 0; k < tpt->size(); ++k)
        if (it < 0 || tpt->at(k) > tpt->at(it)) it = (int) k;
      if (it < 0 || tpt->at(it) < 20 || std::fabs(teta->at(it)) >= 0.7) continue;
      tt->GetEntry(i);
      s.nsel++;

      const double pt = tpt->at(it), eta = teta->at(it), phi = tphi->at(it);
      s.h["vz"]->Fill(vz);
      s.h["tlead_pt"]->Fill(pt);
      s.h["tlead_eta"]->Fill(eta);
      double girth = 0;
      for (size_t c = 0; c < tcpt->at(it).size(); ++c)
      {
        const double dphi = TVector2::Phi_mpi_pi(tcphi->at(it)[c] - phi);
        girth += tcpt->at(it)[c] / pt * std::hypot(tceta->at(it)[c] - eta, dphi);
      }
      s.h["tlead_girth"]->Fill(girth);
      s.h["tlead_nconst"]->Fill(tcpt->at(it).size());

      int ir = -1;
      double best = 0.3;
      for (size_t k = 0; k < rpt->size(); ++k)
      {
        const double dr = std::hypot(reta->at(k) - eta, TVector2::Phi_mpi_pi(rphi->at(k) - phi));
        if (dr < best) { best = dr; ir = (int) k; }
      }
      s.h["match_eff"]->Fill(ir >= 0 ? 1 : 0);
      if (ir >= 0)
      {
        s.h["resp"]->Fill(rpt->at(ir) / pt);
        s.h["deta"]->Fill(reta->at(ir) - eta);
        s.h["rmatch_nconst"]->Fill(rcpt->at(ir).size());
      }
      double rl = 0;
      for (float x : *rpt) rl = std::max(rl, (double) x);
      s.h["rlead_pt"]->Fill(rl);

      double se = 0, si = 0, so = 0;
      int nt = 0;
      for (size_t k = 0; k < em->size(); ++k)
      {
        const double e = std::max(0.f, em->at(k)), a = std::max(0.f, ih->at(k)), b = std::max(0.f, oh->at(k));
        se += e; si += a; so += b;
        if (e + a + b > 0.05) nt++;
        s.h["eprof"]->Fill(etabin->at(k) + 0.5, e + a + b);
      }
      s.h["e_em"]->Fill(se);
      s.h["e_ih"]->Fill(si);
      s.h["e_oh"]->Fill(so);
      s.h["e_tot"]->Fill(se + si + so);
      s.h["ntow"]->Fill(nt);
    }
    delete jt;
    delete tt;
    delete gt;
    for (auto &kv : s.h)
    {
      if (kv.first == "eprof") kv.second->Scale(1.0 / s.nsel);  // mean energy per event per row
      else if (kv.second->Integral() > 0) kv.second->Scale(1.0 / kv.second->Integral());
    }
    printf("%-50s read %7ld, selected %7ld (%.1f%%)%s\n", file.c_str(), s.nread, s.nsel,
           100.0 * s.nsel / s.nread, vzcut > 0 ? Form(" [|vz| < %.0f cm]", vzcut) : "");
    return s;
  }

  void style(TH1D *h, int color, int marker, int lstyle)
  {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.7);
    h->SetLineWidth(2);
    h->SetLineStyle(lstyle);
  }

  void panel(TPad *pad, TH1D *hn, TH1D *ha, TH1D *hc, const char *title, bool logy)
  {
    pad->cd();
    TPad *up = new TPad(Form("u_%s", hn->GetName()), "", 0, 0.32, 1, 1);
    TPad *dn = new TPad(Form("d_%s", hn->GetName()), "", 0, 0, 1, 0.32);
    up->SetBottomMargin(0.02);
    dn->SetTopMargin(0.02);
    dn->SetBottomMargin(0.3);
    up->Draw();
    dn->Draw();
    up->cd();
    up->SetLogy(logy);
    double mx = std::max({hn->GetMaximum(), ha->GetMaximum(), hc->GetMaximum()});
    ha->SetTitle(title);
    ha->GetXaxis()->SetLabelSize(0);
    ha->SetMaximum(logy ? mx * 8 : mx * 1.35);
    if (logy) ha->SetMinimum(std::max(1e-6, 0.5 * std::min(hn->GetMinimum(0), hc->GetMinimum(0))));
    else ha->SetMinimum(0);
    ha->Draw("hist");
    hc->Draw("hist same");
    hn->Draw("e1 same");
    TLegend *l = new TLegend(0.52, 0.70, 0.89, 0.89);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(hn, "new (vertex fixed at 0)", "lp");
    l->AddEntry(ha, "v2 genesis, all vertices", "l");
    l->AddEntry(hc, "v2 genesis, |v_{z}|<10 cm", "l");
    l->Draw();
    dn->cd();
    TH1D *r = (TH1D *) hn->Clone(Form("r_%s", hn->GetName()));
    r->Divide(hc);
    r->SetTitle("");
    r->GetYaxis()->SetTitle("new / v2(|v_{z}|<10)");
    r->GetYaxis()->SetRangeUser(0.5, 1.5);
    r->GetYaxis()->SetNdivisions(505);
    r->GetYaxis()->SetTitleSize(0.09);
    r->GetYaxis()->SetTitleOffset(0.5);
    r->GetYaxis()->SetLabelSize(0.09);
    r->GetXaxis()->SetTitle(hn->GetXaxis()->GetTitle());
    r->GetXaxis()->SetTitleSize(0.11);
    r->GetXaxis()->SetLabelSize(0.09);
    r->Draw("e1");
    TLine *one = new TLine(r->GetXaxis()->GetXmin(), 1, r->GetXaxis()->GetXmax(), 1);
    one->SetLineStyle(2);
    one->Draw();
  }
}  // namespace

void compare_sim(const char *outdir, const char *newVac, const char *newMed, const char *v2Vac,
                 const char *v2Med)
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outdir, true);

  Sample nv = fillSample(newVac, "nv", 0), nm = fillSample(newMed, "nm", 0);
  Sample av = fillSample(v2Vac, "av", 0), am = fillSample(v2Med, "am", 0);
  Sample cv = fillSample(v2Vac, "cv", 10), cm = fillSample(v2Med, "cm", 10);

  TString pdf = Form("%s/compare_sim.pdf", outdir);
  TCanvas *c = new TCanvas("c", "", 1100, 500);
  c->Print(pdf + "[");
  printf("\n%-14s | %-36s | %-36s\n", "", "vacuum: mean new / v2 all / v2 |vz|<10", "medium: mean new / v2 all / v2 |vz|<10");
  int page = 0;
  for (const auto &o : kObs)
  {
    TH1D *a[6] = {nv.h[o.name], av.h[o.name], cv.h[o.name], nm.h[o.name], am.h[o.name], cm.h[o.name]};
    style(a[0], kBlue + 1, 20, 1); style(a[1], kGray + 1, 1, 1); style(a[2], kBlue + 1, 1, 2);
    style(a[3], kRed + 1, 20, 1);  style(a[4], kGray + 1, 1, 1); style(a[5], kRed + 1, 1, 2);
    c->Clear();
    c->Divide(2, 1);
    panel((TPad *) c->cd(1), a[0], a[1], a[2], "vacuum", o.logy);
    panel((TPad *) c->cd(2), a[3], a[4], a[5], "medium", o.logy);
    c->Print(pdf);
    c->Print(Form("%s/page%02d.png", outdir, ++page));
    printf("%-14s | %9.4g / %9.4g / %9.4g     | %9.4g / %9.4g / %9.4g\n", o.name.c_str(),
           a[0]->GetMean(), a[1]->GetMean(), a[2]->GetMean(), a[3]->GetMean(), a[4]->GetMean(), a[5]->GetMean());
  }
  c->Print(pdf + "]");
  printf("\nwrote %s\n", pdf.Data());
}
