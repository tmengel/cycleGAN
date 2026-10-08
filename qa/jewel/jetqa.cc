// jetqa.cc
//
// Jet-level QA histograms from a JEWEL HepMC2 (IO_GenEvent) file, for
// comparing productions. Final-state particles (status 1, excluding
// |pdg| 12-16, as in jetFilter / TruthJetInput) are clustered anti-kT R=0.4.
// An event is SELECTED if its leading jet within |eta| < 0.7 has pT > 20 GeV,
// which is exactly the genesis filter, so filtered and unfiltered samples can
// be compared on equal footing.
//
// Every histogram is filled twice: "w_" with the JEWEL event weight (physical
// distributions) and "u_" unweighted (what a training sample sees). The
// counters in "counts" are: [0] events read, [1] events selected,
// [2] sum of weights read, [3] sum of weights selected.
//
// Usage: jetqa <input.hepmc> <output.root> [max_events (-1 = all)] [leading-jet ptcut, default 20; 0 = no selection]
// Build: see validation/run_compare.sh

#include "fastjet/ClusterSequence.hh"

#include <TFile.h>
#include <TH1D.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const double R = 0.4, ETAMAX = 0.7;
double PTCUT = 20.0;  // leading-jet selection; argument 4 overrides (0 = no selection)

struct HPair {
    TH1D* w;
    TH1D* u;
    void fill(double x, double weight, double extra = 1.0) {
        w->Fill(x, weight * extra);
        u->Fill(x, extra);
    }
};

std::map<std::string, HPair> H;

void book(const std::string& name, const std::string& title, int n, double lo, double hi) {
    H[name] = {new TH1D(("w_" + name).c_str(), title.c_str(), n, lo, hi),
               new TH1D(("u_" + name).c_str(), title.c_str(), n, lo, hi)};
    H[name].w->Sumw2();
    H[name].u->Sumw2();
}

double eventWeight(const std::string& line) {
    std::istringstream iss(line);
    std::string tag;
    long ival;
    double dval;
    iss >> tag >> ival >> ival >> dval >> dval >> dval >> ival >> ival >> ival >> ival >> ival;
    int n = 0;
    iss >> n;
    for (int i = 0; i < n; ++i) iss >> ival;
    iss >> n;
    double w = 1.0;
    if (n > 0) iss >> w;
    return w;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <input.hepmc> <output.root> [max_events]\n";
        return 1;
    }
    const long maxEvents = argc > 3 ? std::atol(argv[3]) : -1;
    if (argc > 4) PTCUT = std::atof(argv[4]);

    TFile out(argv[2], "RECREATE");
    book("lead_pt", ";leading jet p_{T} [GeV]", 40, 0, 80);
    book("lead_eta", ";leading jet #eta", 28, -0.7, 0.7);
    book("lead_nconst", ";leading jet constituents", 50, 0, 50);
    book("lead_mass", ";leading jet mass [GeV]", 40, 0, 20);
    book("lead_girth", ";leading jet girth g", 40, 0, 0.2);
    book("lead_z", ";constituent z = p_{T}/p_{T}^{jet} (leading jet)", 50, 0, 1);
    book("lead_profile", ";r (constituent distance to leading jet axis);#sum p_{T}(r)/p_{T}^{jet}", 8, 0, 0.4);
    book("sublead_pt", ";subleading jet p_{T} [GeV] (|#eta|<0.7)", 35, 0, 70);
    book("njets10", ";jets with p_{T}>10 GeV, |#eta|<0.7", 6, 0, 6);
    book("nfinal", ";final-state particles per event", 50, 0, 200);
    book("log10w", ";log_{10}(event weight)", 60, -16, -4);
    TH1D* counts = new TH1D("counts", "read, selected, sumw read, sumw selected", 4, 0, 4);

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence::set_fastjet_banner_stream(nullptr);

    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open " << argv[1] << "\n";
        return 1;
    }

    std::vector<fastjet::PseudoJet> parts;
    double weight = 1.0;
    long nRead = 0, nSel = 0;
    double swRead = 0, swSel = 0;
    bool inEvent = false;

    auto process = [&]() {
        if (!inEvent) return;
        nRead++;
        swRead += weight;
        fastjet::ClusterSequence cs(parts, jetDef);
        std::vector<fastjet::PseudoJet> jets;
        for (auto& j : fastjet::sorted_by_pt(cs.inclusive_jets(5.0)))
            if (std::fabs(j.eta()) < ETAMAX) jets.push_back(j);
        if (jets.empty() || jets[0].pt() < PTCUT) return;
        nSel++;
        swSel += weight;

        const fastjet::PseudoJet& lj = jets[0];
        H["lead_pt"].fill(lj.pt(), weight);
        H["lead_eta"].fill(lj.eta(), weight);
        H["lead_mass"].fill(lj.m(), weight);
        H["sublead_pt"].fill(jets.size() > 1 ? jets[1].pt() : 0.0, weight);
        int n10 = 0;
        for (auto& j : jets) if (j.pt() > 10) n10++;
        H["njets10"].fill(n10, weight);
        H["nfinal"].fill(parts.size(), weight);
        H["log10w"].fill(std::log10(weight), weight);

        auto constituents = lj.constituents();
        H["lead_nconst"].fill(constituents.size(), weight);
        double girth = 0;
        for (auto& c : constituents) {
            double dr = c.delta_R(lj), z = c.pt() / lj.pt();
            girth += z * dr;
            H["lead_z"].fill(z, weight);
            H["lead_profile"].fill(dr, weight, z);
        }
        H["lead_girth"].fill(girth, weight);
    };

    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("E ", 0) == 0) {
            process();
            if (maxEvents > 0 && nRead >= maxEvents) { inEvent = false; break; }
            parts.clear();
            weight = eventWeight(line);
            inEvent = true;
            continue;
        }
        if (!inEvent || line.rfind("P ", 0) != 0) continue;
        std::istringstream iss(line);
        std::string t;
        long bc;
        int id, st;
        double px, py, pz, e, m;
        iss >> t >> bc >> id >> px >> py >> pz >> e >> m >> st;
        if (st == 1 && !(std::abs(id) >= 12 && std::abs(id) <= 16)) parts.emplace_back(px, py, pz, e);
    }
    process();

    counts->SetBinContent(1, nRead);
    counts->SetBinContent(2, nSel);
    counts->SetBinContent(3, swRead);
    counts->SetBinContent(4, swSel);
    out.Write();
    std::cout << argv[1] << ": read " << nRead << ", selected " << nSel
              << " (" << (nRead ? 100.0 * nSel / nRead : 0) << "%), sumw read " << swRead
              << ", sumw selected " << swSel << "\n";
    return 0;
}
