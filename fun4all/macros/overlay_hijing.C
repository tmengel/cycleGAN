// overlay_hijing.C -- add a HIJING 0-10% background image to every JEWEL image.
//
//   root -b -q 'overlay_hijing.C("<jobs list>", "<jewel images dir>", "<out dir>",
//                                "<sample>", <run>, "<hijing dir>", "<hijing counts>")'
//
// <jobs list>: lines "<segment> <nevents> <first global event>" (written by
// scripts/submit_overlay.sh from the pass 1 job list). For each segment it reads
// <jewel images dir>/images_<sample>-<run>-<segment>.root (histograms
// h_eta_phi_cent0_file<seg>_evt<i>, i = 0..nevents-1) and writes
// <out dir>/images_<sample>_hijing-<run>-<segment>.root with the same histogram names, where
//   overlay(event i) = JEWEL image(event i) + HIJING image(global event e = first + i).
// Image contents are non-negative tower energies on the same 24x64 grid in both inputs, so the
// sum equals yeonjugo's combine_pythia_hijing_v3_jewel.C (tower-level sum) followed by imaging.
//
// HIJING event e = image number p (key order) of file N (see hijing/README.md). Events with
// e >= number of HIJING images are not written (no background reuse); a segment that ends up
// with no images is not written at all.
//
// Each output file also holds a TTree "hijing_map" (evt, global_event, hijing_file,
// hijing_key) recording which background went into each image.

#include <TFile.h>
#include <TH2D.h>
#include <TKey.h>
#include <TList.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{
  struct HijingIndex
  {
    std::vector<int> file;   // file number of each HIJING file with images, in order
    std::vector<long> first; // global index of its first image
    long total{0};
  };

  HijingIndex readCounts(const std::string &countsFile)
  {
    HijingIndex idx;
    std::ifstream in(countsFile);
    int n, c;
    while (in >> n >> c)
    {
      if (c <= 0) continue;
      idx.file.push_back(n);
      idx.first.push_back(idx.total);
      idx.total += c;
    }
    return idx;
  }

  // keeps the currently needed HIJING file open
  struct HijingReader
  {
    std::string dir;
    TFile *f{nullptr};
    int fileNumber{-1};
    std::vector<std::string> keys;

    TH2 *get(int n, long p, std::string &keyName)
    {
      if (n != fileNumber)
      {
        if (f) f->Close();
        delete f;
        f = TFile::Open(Form("%s/type4_run19_hijingAll_noNoise_images_cent0_file%d.root", dir.c_str(), n));
        if (!f || f->IsZombie())
        {
          std::cerr << "overlay_hijing: cannot open HIJING file " << n << std::endl;
          return nullptr;
        }
        fileNumber = n;
        keys.clear();
        for (TObject *o : *f->GetListOfKeys()) keys.push_back(o->GetName());
      }
      if (p < 0 || p >= (long) keys.size()) return nullptr;
      keyName = keys[p];
      return dynamic_cast<TH2 *>(f->Get(keyName.c_str()));
    }
  };
}  // namespace

int overlay_hijing(const std::string &jobsList, const std::string &jewelDir, const std::string &outDir,
                   const std::string &sample, const int run, const std::string &hijingDir,
                   const std::string &countsFile)
{
  const HijingIndex idx = readCounts(countsFile);
  if (idx.total == 0)
  {
    std::cerr << "overlay_hijing: no HIJING images in " << countsFile << std::endl;
    return 1;
  }
  HijingReader hij;
  hij.dir = hijingDir;

  std::ifstream jobs(jobsList);
  long seg, nev, firstEvent;
  int nfiles = 0;
  long nimages = 0;
  while (jobs >> seg >> nev >> firstEvent)
  {
    if (firstEvent >= idx.total) continue;
    TString jewelName = Form("%s/images_%s-%010d-%06ld.root", jewelDir.c_str(), sample.c_str(), run, seg);
    TFile *jf = TFile::Open(jewelName);
    if (!jf || jf->IsZombie())
    {
      std::cerr << "overlay_hijing: cannot open " << jewelName << std::endl;
      return 1;
    }
    TString outName = Form("%s/images_%s_hijing-%010d-%06ld.root", outDir.c_str(), sample.c_str(), run, seg);
    TFile *out = new TFile(outName, "RECREATE");
    int evt = 0;
    long gevt = 0;
    int hfile = 0;
    char hkey[256];
    TTree *map = new TTree("hijing_map", "HIJING background used for each image");
    map->Branch("evt", &evt);
    map->Branch("global_event", &gevt);
    map->Branch("hijing_file", &hfile);
    map->Branch("hijing_key", hkey, "hijing_key/C");

    int written = 0;
    for (int i = 0; i < nev; ++i)
    {
      gevt = firstEvent + i;
      if (gevt >= idx.total) break;  // no reuse
      TString name = Form("h_eta_phi_cent0_file%ld_evt%d", seg, i);
      TH2 *hj = dynamic_cast<TH2 *>(jf->Get(name));
      if (!hj)
      {
        std::cerr << "overlay_hijing: " << name << " missing in " << jewelName << std::endl;
        return 1;
      }
      // locate HIJING event gevt: last file whose first index <= gevt
      const size_t k = std::upper_bound(idx.first.begin(), idx.first.end(), gevt) - idx.first.begin() - 1;
      std::string key;
      TH2 *hh = hij.get(idx.file[k], gevt - idx.first[k], key);
      if (!hh || hh->GetNbinsX() != hj->GetNbinsX() || hh->GetNbinsY() != hj->GetNbinsY())
      {
        std::cerr << "overlay_hijing: bad HIJING image for global event " << gevt << std::endl;
        return 1;
      }
      // Bin-by-bin sum over the 24x64 tower grid. TH2::Add refuses these histograms: the HIJING
      // images carry draw_event_v2.C's phi axis [0, 2*3.14156], ours the exact [0, 2 pi]; the
      // bins are the same towers. The output keeps the JEWEL image's axes.
      TH2 *sum = (TH2 *) hj->Clone(name);
      sum->SetDirectory(nullptr);
      for (int bx = 1; bx <= hj->GetNbinsX(); ++bx)
        for (int by = 1; by <= hj->GetNbinsY(); ++by)
          sum->SetBinContent(bx, by, hj->GetBinContent(bx, by) + hh->GetBinContent(bx, by));
      out->cd();
      sum->Write();
      delete sum;
      delete hh;
      evt = i;
      hfile = idx.file[k];
      snprintf(hkey, sizeof(hkey), "%s", key.c_str());
      map->Fill();
      ++written;
    }
    out->cd();
    map->Write();
    out->Close();
    jf->Close();
    if (written == 0) gSystem->Unlink(outName);
    else { ++nfiles; nimages += written; }
  }
  std::cout << "overlay_hijing: wrote " << nimages << " images in " << nfiles << " files" << std::endl;
  return 0;
}
