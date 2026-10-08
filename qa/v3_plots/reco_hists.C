// reco_hists.C -- reco-level QA histograms for one slice of a v3 sample's trees (one condor job).
//
//   root -b -q 'reco_hists.C("<sample>", <run>, <first segment>, <last segment>, "<out.root>")'
//
// Same selection and definitions as the slide-6 figures in make_figs.C: events with a leading
// truth jet (anti-kT R = 0.4) pT > 20 GeV, |eta| < 0.7; leading reco jet pT; eta of all reco jets
// (pT >= 10 GeV, the tree cut); reco jet matched to the leading truth jet within dR < 0.3.
// Missing segment files are skipped (reported at the end).

const char *TAG = "/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v3_2026-10-03";

void reco_hists(const char *sample, int run, int first, int last, const char *outfile)
{
  TChain jt("jettree");
  int nmissing = 0;
  for (int i = first; i <= last; ++i)
  {
    TString f = Form("%s/trees/%s/trees_%s-%010d-%06d.root", TAG, sample, sample, run, i);
    if (gSystem->AccessPathName(f)) { nmissing++; continue; }
    jt.Add(f);
  }
  std::vector<float> *tp = nullptr, *te = nullptr, *tf = nullptr, *rp = nullptr, *re = nullptr, *rf = nullptr;
  jt.SetBranchAddress("tjet_pt", &tp);
  jt.SetBranchAddress("tjet_eta", &te);
  jt.SetBranchAddress("tjet_phi", &tf);
  jt.SetBranchAddress("rjet_pt", &rp);
  jt.SetBranchAddress("rjet_eta", &re);
  jt.SetBranchAddress("rjet_phi", &rf);

  TFile out(outfile, "RECREATE");
  TH1D tpt("tpt", ";leading truth jet p_{T} [GeV];events (normalised)", 30, 0, 60);
  TH1D rpt("rpt", ";leading reco jet p_{T} [GeV];events (normalised)", 30, 0, 60);
  TH1D eta("eta", ";reco jet #eta (p_{T} > 10 GeV);jets (normalised)", 24, -1.2, 1.2);
  TH2D resp("resp", ";leading truth jet p_{T} [GeV];matched reco jet p_{T} [GeV]", 30, 0, 60, 30, 0, 60);
  TH1D counts("counts", "events read, events selected", 2, 0, 2);

  const Long64_t n = jt.GetEntries();
  long nsel = 0;
  for (Long64_t i = 0; i < n; ++i)
  {
    jt.GetEntry(i);
    int it = -1;
    for (size_t k = 0; k < tp->size(); ++k) if (it < 0 || tp->at(k) > tp->at(it)) it = k;
    // no event selection: every event is used (jets in the trees are >= 10 GeV)
    nsel++;
    if (it >= 0) tpt.Fill(tp->at(it));
    double rl = 0, best = 0.3;
    int im = -1;
    for (size_t k = 0; k < rp->size(); ++k)
    {
      rl = std::max(rl, (double) rp->at(k));
      eta.Fill(re->at(k));
      if (it < 0) continue;
      const double dr = std::hypot(re->at(k) - te->at(it), TVector2::Phi_mpi_pi(rf->at(k) - tf->at(it)));
      if (dr < best) { best = dr; im = k; }
    }
    if (rl > 0) rpt.Fill(rl);
    if (it >= 0 && im >= 0) resp.Fill(tp->at(it), rp->at(im));
  }
  counts.SetBinContent(1, n);
  counts.SetBinContent(2, nsel);
  out.Write();
  out.Close();
  printf("%s segments %d-%d: %lld events, %ld selected, %d files missing\n", sample, first, last, n, nsel, nmissing);
}
