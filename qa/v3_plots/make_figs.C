// make_figs.C -- figures for cyclegan_update_v3.pptx (v3_2026-10-03 dataset).
//
//   root -b -q 'make_figs.C("figs")'           (or a subset: make_figs.C("figs", "gen,evt,reco"))
//
// gen_*.png   generator level, vacuum vs no-recoil: data/gen_<sample>.root = qa/jewel/jetqa
//             on the tag's hepmc/<sample>_all.hepmc (every simulated event, ~1M per sample);
//             event weights applied; no event selection (leading anti-kT R=0.4 jet in |eta| < 0.7)
// evt_*.png   one calorimeter image from each of the four samples
// reco_*.png  truth vs reco jets from all v3 trees (~1M events per sample): data/reco_<sample>.root,
//             the hadd of the condor outputs of reco_hists.C (data/reco_hists.job); events with a
//             no event selection (jets in the trees are >= 10 GeV)

const char *TAG = "/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v3_2026-10-03";
const char *DATA = "/sphenix/user/tmengel/cycleGAN/qa/v3_plots/data";

void label(const char *line1, const char *line2 = nullptr, double x = 0.21, double y = 0.86)
{
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.045);
  l.DrawLatex(x, y, "#bf{#it{sPHENIX}} Simulation");
  l.SetTextSize(0.038);
  l.DrawLatex(x, y - 0.055, line1);
  if (line2) l.DrawLatex(x, y - 0.105, line2);
}

void styleh(TH1 *h, int color, int ls)
{
  h->SetLineColor(color);
  h->SetLineWidth(3);
  h->SetLineStyle(ls);
  h->SetStats(0);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetTitleOffset(1.3);
}

void gen_figs(const char *out)
{
  TFile *fv = TFile::Open(Form("%s/gen_vacuum.root", DATA));
  TFile *fm = TFile::Open(Form("%s/gen_medium.root", DATA));
  struct P { const char *h, *x; bool logy; };
  std::vector<P> ps = {{"w_lead_pt", "leading jet p_{T} [GeV]", true},
                       {"w_lead_eta", "leading jet #eta", false},
                       {"w_lead_mass", "leading jet mass [GeV]", false},
                       {"w_lead_girth", "leading jet girth", false},
                       {"w_lead_nconst", "leading jet constituents", false},
                       {"w_nfinal", "final-state particles / event", false}};
  for (auto &p : ps)
  {
    TCanvas c("c", "", 700, 560);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetLogy(p.logy);
    TH1 *v = (TH1 *) fv->Get(p.h)->Clone("v"), *m = (TH1 *) fm->Get(p.h)->Clone("m");
    v->Scale(1.0 / v->Integral());
    m->Scale(1.0 / m->Integral());
    styleh(v, kAzure + 2, 1);
    styleh(m, kOrange + 7, 2);
    v->SetTitle(Form(";%s;normalised", p.x));
    double mx = std::max(v->GetMaximum(), m->GetMaximum());
    v->SetMaximum(p.logy ? mx * 3000 : mx * 1.55);
    if (!p.logy) v->SetMinimum(0);
    v->Draw("hist");
    m->Draw("hist same");
    label("JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV", "generator level, weighted");
    TLegend lg(0.62, 0.72, 0.9, 0.88);
    lg.SetBorderSize(0);
    lg.SetFillStyle(0);
    lg.SetTextSize(0.04);
    const bool showMean = TString(p.h) != "w_lead_eta";
    lg.AddEntry(v, showMean ? Form("vacuum, mean %.3g", v->GetMean()) : "vacuum", "l");
    lg.AddEntry(m, showMean ? Form("no-recoil, mean %.3g", m->GetMean()) : "no-recoil", "l");
    lg.Draw();
    c.SaveAs(Form("%s/gen_%s.png", out, p.h + 2));
  }
}

void evt_figs(const char *out)
{
  // one image per sample; vacuum and vacuum+HIJING are the same JEWEL event (likewise medium)
  struct S { const char *dir, *file, *title; int evt; };
  std::vector<S> ss = {{"vacuum", "images_vacuum-0000001000-000100.root", "JEWEL vacuum", 1},
                       {"medium", "images_medium-0000002000-000100.root", "JEWEL no-recoil", 1},
                       {"vacuum_hijing", "images_vacuum_hijing-0000001000-000100.root", "JEWEL vacuum + HIJING 0-10%", 1},
                       {"medium_hijing", "images_medium_hijing-0000002000-000100.root", "JEWEL no-recoil + HIJING 0-10%", 1}};
  gStyle->SetPalette(kBird);
  for (int j : {0, 1})
  {
    TFile *f = TFile::Open(Form("%s/images/%s/%s", TAG, ss[j].dir, ss[j].file));
    int best = 0;
    double emax = -1;
    for (int e = 0; e < 10; ++e)
    {
      TH2 *h = (TH2 *) f->Get(Form("h_eta_phi_cent0_file100_evt%d", e));
      if (h && h->Integral() > emax) { emax = h->Integral(); best = e; }
    }
    ss[j].evt = best;
    ss[j + 2].evt = best;
    f->Close();
  }
  for (auto &s : ss)
  {
    TFile *f = TFile::Open(Form("%s/images/%s/%s", TAG, s.dir, s.file));
    if (!f || f->IsZombie()) { printf("missing %s\n", s.file); continue; }
    TH2 *h = (TH2 *) f->Get(Form("h_eta_phi_cent0_file100_evt%d", s.evt));
    TCanvas c("c", "", 700, 620);
    c.SetRightMargin(0.17);
    c.SetLeftMargin(0.13);
    c.SetBottomMargin(0.13);
    h->SetStats(0);
    h->SetTitle(";#eta;#phi [rad]");
    h->GetZaxis()->SetTitle("E [GeV]");
    h->GetZaxis()->SetTitleOffset(1.2);
    h->SetMinimum(1e-3);
    h->Draw("colz");
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.042);
    l.DrawLatex(0.15, 0.93, Form("%s  (#SigmaE = %.0f GeV)", s.title, h->Integral()));
    c.SaveAs(Form("%s/evt_%s.png", out, s.dir));
  }
}

void reco_figs(const char *out)
{
  struct S { const char *name, *title; int run, color, ls; };
  std::vector<S> ss = {{"vacuum", "JEWEL vacuum", 1000, kAzure + 2, 1}, {"medium", "JEWEL no-recoil", 2000, kOrange + 7, 2}};
  std::map<std::string, TH1D *> tpt, rpt, eta;
  std::map<std::string, TH2D *> resp;
  for (auto &s : ss)
  {
    std::string n = s.name;
    TFile *f = TFile::Open(Form("%s/reco_%s.root", DATA, s.name));
    tpt[n] = (TH1D *) f->Get("tpt");
    rpt[n] = (TH1D *) f->Get("rpt");
    eta[n] = (TH1D *) f->Get("eta");
    resp[n] = (TH2D *) f->Get("resp");
    TH1D *c = (TH1D *) f->Get("counts");
    printf("%s: %.0f events read, %.0f selected\n", s.name, c->GetBinContent(1), c->GetBinContent(2));
  }
  auto pair = [&](std::map<std::string, TH1D *> &hm, const char *fname, bool logy) {
    TCanvas c("c", "", 700, 560);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetLogy(logy);
    TH1D *v = hm["vacuum"], *m = hm["medium"];
    v->Scale(1.0 / v->Integral());
    m->Scale(1.0 / m->Integral());
    styleh(v, kAzure + 2, 1);
    styleh(m, kOrange + 7, 2);
    double mx = std::max(v->GetMaximum(), m->GetMaximum());
    v->SetMaximum(logy ? mx * 300 : mx * 1.55);
    if (logy) v->SetMinimum(1e-4); else v->SetMinimum(0);
    v->Draw("hist");
    m->Draw("hist same");
    label("JEWEL #sqrt{s_{NN}} = 200 GeV, full GEANT4", "all events, no selection");
    TLegend lg(0.66, 0.58, 0.9, 0.70);
    lg.SetBorderSize(0);
    lg.SetFillStyle(0);
    lg.SetTextSize(0.04);
    lg.AddEntry(v, "vacuum", "l");
    lg.AddEntry(m, "no-recoil", "l");
    lg.Draw();
    c.SaveAs(Form("%s/%s.png", out, fname));
  };
  pair(tpt, "reco_truth_pt", true);
  pair(rpt, "reco_reco_pt", true);
  pair(eta, "reco_eta", false);
  gStyle->SetPalette(kBird);
  for (auto &s : ss)
  {
    TCanvas c("c", "", 700, 620);
    c.SetRightMargin(0.15);
    c.SetLeftMargin(0.14);
    c.SetBottomMargin(0.13);
    TH2D *h = resp[s.name];
    c.SetLogz();
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->Draw("colz");
    TLine d(0, 0, 60, 60);
    d.SetLineStyle(2);
    d.Draw();
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.17, 0.86, s.title);
    c.SaveAs(Form("%s/reco_resp_%s.png", out, s.name));
  }
}

void make_figs(const char *out = "figs", const char *which = "gen,evt,reco")
{
  gSystem->mkdir(out, true);
  gROOT->SetBatch(kTRUE);
  TString w = which;
  if (w.Contains("gen")) gen_figs(out);
  if (w.Contains("evt")) evt_figs(out);
  if (w.Contains("reco")) reco_figs(out);
}
