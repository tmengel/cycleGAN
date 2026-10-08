// compare.C -- overlay jetqa outputs of the new production against Luis's
// genesis sample, for vacuum and medium.
//
// root -b -q 'compare.C("compare_luis", 100000, 100000, 194042, 271545)'
//   ngen_*: events GENERATED (before any jet filter) behind each sample; used
//   to normalise yields per generated event for the medium/vacuum ratio.
//
// Writes <dir>/compare.pdf and prints a summary table (means, shape chi2).

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>

#include <cstdio>
#include <vector>

namespace {

TH1D* get(TFile* f, const char* name, const char* tag) {
  TH1D* h = (TH1D*) f->Get(name)->Clone(Form("%s_%s", name, tag));
  h->SetDirectory(nullptr);
  return h;
}

void style(TH1D* h, int color, int marker) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.7);
  h->SetLineWidth(2);
}

// Upper pad: new (points) vs Luis (line); lower pad: new / Luis.
void overlay(TPad* pad, TH1D* hn, TH1D* hl, const char* title, bool logy) {
  pad->cd();
  TPad* up = new TPad(Form("up_%s", hn->GetName()), "", 0, 0.32, 1, 1);
  TPad* dn = new TPad(Form("dn_%s", hn->GetName()), "", 0, 0, 1, 0.32);
  up->SetBottomMargin(0.02);
  dn->SetTopMargin(0.02);
  dn->SetBottomMargin(0.3);
  up->Draw();
  dn->Draw();

  up->cd();
  up->SetLogy(logy);
  hl->SetTitle(title);
  hl->GetXaxis()->SetLabelSize(0);
  double mx = std::max(hn->GetMaximum(), hl->GetMaximum());
  hl->SetMaximum(logy ? mx * 5 : mx * 1.3);
  if (logy) hl->SetMinimum(std::max(1e-6 * mx, 0.3 * hl->GetMinimum(0)));
  hl->Draw("hist");
  hn->Draw("e1 same");
  TLegend* leg = new TLegend(0.55, 0.72, 0.89, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(hn, "new (this setup)", "lp");
  leg->AddEntry(hl, "Luis genesis", "l");
  leg->Draw();

  dn->cd();
  TH1D* r = (TH1D*) hn->Clone(Form("%s_ratio", hn->GetName()));
  r->Divide(hl);
  r->SetTitle("");
  r->GetYaxis()->SetTitle("new / Luis");
  r->GetYaxis()->SetRangeUser(0.5, 1.5);
  r->GetYaxis()->SetNdivisions(505);
  r->GetYaxis()->SetTitleSize(0.1);
  r->GetYaxis()->SetTitleOffset(0.45);
  r->GetYaxis()->SetLabelSize(0.09);
  r->GetXaxis()->SetTitle(hn->GetXaxis()->GetTitle());
  r->GetXaxis()->SetTitleSize(0.11);
  r->GetXaxis()->SetLabelSize(0.09);
  r->Draw("e1");
  TLine* one = new TLine(r->GetXaxis()->GetXmin(), 1, r->GetXaxis()->GetXmax(), 1);
  one->SetLineStyle(2);
  one->Draw();
}

double chi2ndf(TH1D* a, TH1D* b) {
  // shape test on weighted histograms
  return a->Chi2Test(b, "WW CHI2/NDF");
}

}  // namespace

void compare(const char* dir, double ngenNewVac, double ngenNewMed, double ngenLuisVac,
             double ngenLuisMed) {
  gStyle->SetOptStat(0);
  TFile* fnv = TFile::Open(Form("%s/new_vacuum.root", dir));
  TFile* fnm = TFile::Open(Form("%s/new_medium.root", dir));
  TFile* flv = TFile::Open(Form("%s/luis_vacuum.root", dir));
  TFile* flm = TFile::Open(Form("%s/luis_norecoil.root", dir));
  TString pdf = Form("%s/compare.pdf", dir);

  const std::vector<const char*> obs = {"lead_pt",     "lead_eta",   "lead_nconst", "lead_mass",
                                        "lead_girth",  "lead_z",     "lead_profile", "sublead_pt",
                                        "njets10",     "nfinal",     "log10w"};
  const std::vector<bool> logy = {true, false, false, true, false, true, true, true, false, false, true};

  TCanvas* c = new TCanvas("c", "", 1100, 500);
  int page = 0;
  c->Print(pdf + "[");

  // ---- page 1: absolute yields per generated event, + medium/vacuum ratio
  {
    TH1D* nv = get(fnv, "w_lead_pt", "nv");
    TH1D* nm = get(fnm, "w_lead_pt", "nm");
    TH1D* lv = get(flv, "w_lead_pt", "lv");
    TH1D* lm = get(flm, "w_lead_pt", "lm");
    nv->Scale(1.0 / ngenNewVac, "width");
    nm->Scale(1.0 / ngenNewMed, "width");
    lv->Scale(1.0 / ngenLuisVac, "width");
    lm->Scale(1.0 / ngenLuisMed, "width");
    for (auto* h : {nv, nm, lv, lm}) h->GetYaxis()->SetTitle("#sum w / N_{gen} / GeV");
    style(nv, kBlue + 1, 20); style(lv, kBlue + 1, 1);
    style(nm, kRed + 1, 20);  style(lm, kRed + 1, 1);
    c->Clear();
    c->Divide(2, 1);
    overlay((TPad*) c->cd(1), nv, lv, "vacuum: weighted leading-jet yield per generated event", true);
    overlay((TPad*) c->cd(2), nm, lm, "medium: weighted leading-jet yield per generated event", true);
    c->Print(pdf);
    c->Print(Form("%s/page%02d.png", dir, ++page));

    TH1D* rn = (TH1D*) nm->Clone("raa_new");
    rn->Divide(nv);
    TH1D* rl = (TH1D*) lm->Clone("raa_luis");
    rl->Divide(lv);
    style(rn, kBlack, 20);
    style(rl, kMagenta + 1, 1);
    c->Clear();
    c->cd();
    rn->SetTitle("medium / vacuum leading-jet yield (R_{AA}-like, same N_{gen} normalisation);leading jet p_{T} [GeV];medium / vacuum");
    rn->GetYaxis()->SetRangeUser(0, 1.2);
    rn->Draw("e1");
    rl->Draw("e1 same");
    TLegend* leg = new TLegend(0.6, 0.75, 0.89, 0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(rn, "new (ETAMAX 2.0, own tables)", "lp");
    leg->AddEntry(rl, "Luis genesis", "lp");
    leg->Draw();
    c->Print(pdf);
    c->Print(Form("%s/page%02d.png", dir, ++page));

    printf("\nleading-jet yield ratio medium/vacuum (sum over 20-80 GeV): new %.3f  Luis %.3f\n",
           nm->Integral("width") / nv->Integral("width"), lm->Integral("width") / lv->Integral("width"));
    printf("%-14s %10s %10s %10s %10s\n", "pT bin", "new", "+-", "Luis", "+-");
    for (int b = 1; b <= rn->GetNbinsX(); b += 3) {
      // merge 3 bins for a readable table
      double n1 = 0, n2 = 0, l1 = 0, l2 = 0, en1 = 0, en2 = 0, el1 = 0, el2 = 0;
      for (int k = b; k < b + 3 && k <= rn->GetNbinsX(); ++k) {
        n1 += nm->GetBinContent(k); n2 += nv->GetBinContent(k);
        l1 += lm->GetBinContent(k); l2 += lv->GetBinContent(k);
        en1 += pow(nm->GetBinError(k), 2); en2 += pow(nv->GetBinError(k), 2);
        el1 += pow(lm->GetBinError(k), 2); el2 += pow(lv->GetBinError(k), 2);
      }
      if (n2 <= 0 || l2 <= 0 || n1 <= 0 || l1 <= 0) continue;
      double rN = n1 / n2, rL = l1 / l2;
      printf("%5.0f-%-5.0f GeV %10.3f %10.3f %10.3f %10.3f\n", nm->GetXaxis()->GetBinLowEdge(b),
             nm->GetXaxis()->GetBinUpEdge(std::min(b + 2, rn->GetNbinsX())), rN,
             rN * sqrt(en1 / (n1 * n1) + en2 / (n2 * n2)), rL, rL * sqrt(el1 / (l1 * l1) + el2 / (l2 * l2)));
    }
  }

  // ---- shape pages (weighted, unit area), vacuum left / medium right
  printf("\n%-13s | %-31s | %-31s\n", "", "vacuum: mean new / Luis, chi2/ndf", "medium: mean new / Luis, chi2/ndf");
  for (size_t i = 0; i < obs.size(); ++i) {
    for (const char* kind : {"w", "u"}) {
      if (TString(kind) == "u" && TString(obs[i]) != "lead_pt") continue;  // unweighted: only pT
      TString name = Form("%s_%s", kind, obs[i]);
      TH1D* nv = get(fnv, name, "nv");
      TH1D* nm = get(fnm, name, "nm");
      TH1D* lv = get(flv, name, "lv");
      TH1D* lm = get(flm, name, "lm");
      double cv = chi2ndf(nv, lv), cm = chi2ndf(nm, lm);
      if (TString(obs[i]) == "lead_profile") {
        // pT fraction per annulus per selected jet: normalise by selected sum of weights
        int bin = TString(kind) == "w" ? 4 : 2;
        nv->Scale(1.0 / ((TH1D*) fnv->Get("counts"))->GetBinContent(bin));
        nm->Scale(1.0 / ((TH1D*) fnm->Get("counts"))->GetBinContent(bin));
        lv->Scale(1.0 / ((TH1D*) flv->Get("counts"))->GetBinContent(bin));
        lm->Scale(1.0 / ((TH1D*) flm->Get("counts"))->GetBinContent(bin));
      } else {
        for (auto* h : {nv, nm, lv, lm}) h->Scale(1.0 / h->Integral());
        for (auto* h : {nv, nm, lv, lm}) h->GetYaxis()->SetTitle("normalised");
      }
      style(nv, kBlue + 1, 20); style(lv, kBlue + 1, 1);
      style(nm, kRed + 1, 20);  style(lm, kRed + 1, 1);
      TString tag = TString(kind) == "w" ? "weighted" : "UNWEIGHTED (raw counts)";
      c->Clear();
      c->Divide(2, 1);
      overlay((TPad*) c->cd(1), nv, lv, Form("vacuum, %s", tag.Data()), logy[i]);
      overlay((TPad*) c->cd(2), nm, lm, Form("medium, %s", tag.Data()), logy[i]);
      c->Print(pdf);
    c->Print(Form("%s/page%02d.png", dir, ++page));
      printf("%-13s | %8.4g / %-8.4g  %6.2f      | %8.4g / %-8.4g  %6.2f\n",
             (TString(kind) + "_" + obs[i]).Data(), nv->GetMean(), lv->GetMean(), cv, nm->GetMean(),
             lm->GetMean(), cm);
    }
  }
  c->Print(pdf + "]");
  printf("\nwrote %s\n", pdf.Data());
}
