// make_images.C -- one calorimeter image (TH2D) per event from a pass 2 tree file.
//
//   root -b -q 'make_images.C("trees.root", "images.root", <filenumber>)'
//
// Same images as yeonjugo's create_calo_images/macro/draw_event_v2.C, the macro used for the
// earlier productions, for the pp-like (no centrality) case:
//   * 24 x 64 histogram, eta in [-1.1, 1.1], phi in [0, 2 pi], filled by tower index
//     (etabin + 1, phibin + 1) from towertree;
//   * content = EMCal (retowered) + iHCal + oHCal energy in GeV, each clipped at 0 from below;
//   * histogram name h_eta_phi_cent0_file<filenumber>_evt<entry>, one per towertree entry, in
//     event order (entry i of the image file = entry i of the trees' jettree/globaltree).
// Differences: input/output paths are arguments instead of being assembled from directory
// names, and the phi axis uses TMath::Pi() (draw_event_v2.C has pi = 3.14156, which shifts
// the axis labels slightly; the bin contents are unaffected).

#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

int make_images(const std::string &treeFile, const std::string &imageFile, const int filenumber)
{
  TFile *in = TFile::Open(treeFile.c_str());
  if (!in || in->IsZombie())
  {
    std::cerr << "make_images: cannot open " << treeFile << std::endl;
    return 1;
  }
  TTree *towertree = dynamic_cast<TTree *>(in->Get("towertree"));
  if (!towertree)
  {
    std::cerr << "make_images: no towertree in " << treeFile << std::endl;
    return 1;
  }

  std::vector<float> *em = nullptr;
  std::vector<float> *ih = nullptr;
  std::vector<float> *oh = nullptr;
  std::vector<int> *phibin = nullptr;
  std::vector<int> *etabin = nullptr;
  towertree->SetBranchAddress("emcalenergy", &em);
  towertree->SetBranchAddress("ihcalenergy", &ih);
  towertree->SetBranchAddress("ohcalenergy", &oh);
  towertree->SetBranchAddress("phibin", &phibin);
  towertree->SetBranchAddress("etabin", &etabin);

  TFile *out = new TFile(imageFile.c_str(), "RECREATE");
  TH2D *h = new TH2D("h_eta_phi", ";#eta;#phi", 24, -1.1, 1.1, 64, 0, 2 * TMath::Pi());
  h->SetDirectory(nullptr);

  const Long64_t nentries = towertree->GetEntries();
  for (Long64_t q = 0; q < nentries; ++q)
  {
    towertree->GetEntry(q);
    h->Reset();
    const size_t ncal = em->size();
    for (size_t i = 0; i < ncal; ++i)
    {
      const double e = std::max(0.f, em->at(i)) + std::max(0.f, ih->at(i)) + std::max(0.f, oh->at(i));
      h->SetBinContent(etabin->at(i) + 1, phibin->at(i) + 1, e);
    }
    h->SetName(Form("h_eta_phi_cent0_file%d_evt%lld", filenumber, q));
    out->cd();
    h->Write();
  }
  out->Close();
  in->Close();
  std::cout << "make_images: wrote " << nentries << " images to " << imageFile << std::endl;
  return 0;
}
