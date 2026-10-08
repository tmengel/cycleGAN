// jetFilter.cc
//
// Optional per-event filter for JEWEL HepMC2 (IO_GenEvent) output. Adapted
// from Luis's genesisJetFilter.cc (reference_luis/), same selection:
//
//   final-state particles (status==1), excluding |pdgid| in [12,16]
//   (neutrinos, mu, tau -- matches coresoftware TruthJetInput), clustered
//   anti-kT with radius R. An event PASSES if any jet has pT > ptmin and
//   |eta| < etamax.
//
// Differences from genesisJetFilter:
//   * the output is a complete, standalone HepMC2 file (header + footer), so
//     every filtered chunk can be read on its own or simply concatenated with
//     scripts/merge_hepmc.sh;
//   * an input without the END_EVENT_LISTING footer (truncated JEWEL run) is
//     rejected here instead of relying on the caller to check.
//
// Passing events are copied verbatim (every line from "E " up to the next
// event), so weights, vertices, and cross-section lines are kept untouched.
//
// A sidecar "<output>.stats" file records TOTAL/PASSED counts and the sum of
// event weights before/after the cut (JEWEL events are weighted by default;
// the weight is the first entry of the "E " line's weight list).
//
// Usage:
//   jetFilter <input.hepmc> <output.hepmc> <R> <ptmin_GeV> <etamax>

#include "fastjet/ClusterSequence.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool isExcludedLepton(int pdg) {
    int apdg = std::abs(pdg);
    return apdg >= 12 && apdg <= 16;
}

bool startsWith(const std::string& line, const char* tag) {
    return line.rfind(tag, 0) == 0;
}

// HepMC2 "E" line: E evtno nMPI scale aQCD aQED sigvtx bcvtx nvtx bc1 bc2
//                    nrandom [rnd...] nweights [w...]
double eventWeight(const std::string& line) {
    std::istringstream iss(line);
    std::string tag;
    long ival;
    double dval;
    iss >> tag >> ival >> ival >> dval >> dval >> dval >> ival >> ival >> ival >> ival >> ival;
    int nrnd = 0;
    iss >> nrnd;
    for (int i = 0; i < nrnd; ++i) iss >> ival;
    int nw = 0;
    iss >> nw;
    double w = 1.0;
    if (nw > 0 && (iss >> w)) return w;
    return 1.0;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.hepmc> <output.hepmc> <R> <ptmin_GeV> <etamax>\n";
        return 1;
    }
    const std::string inPath  = argv[1];
    const std::string outPath = argv[2];
    const double R      = std::stod(argv[3]);
    const double ptmin  = std::stod(argv[4]);
    const double etamax = std::stod(argv[5]);

    std::ifstream in(inPath);
    if (!in) {
        std::cerr << "ERROR: cannot open input file " << inPath << "\n";
        return 1;
    }
    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "ERROR: cannot open output file " << outPath << "\n";
        return 1;
    }
    out << "\nHepMC::Version 2.06.05\n"
        << "HepMC::IO_GenEvent-START_EVENT_LISTING\n";

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);

    long nTotal = 0, nPassed = 0;
    double sumwTotal = 0, sumwPassed = 0;
    bool inEvent = false, sawFooter = false;
    double weight = 1.0;
    std::vector<std::string> buf;
    std::vector<fastjet::PseudoJet> parts;

    auto flushEvent = [&]() {
        if (!inEvent) return;
        nTotal++;
        sumwTotal += weight;
        bool pass = false;
        if (!parts.empty()) {
            fastjet::ClusterSequence cs(parts, jetDef);
            for (const fastjet::PseudoJet& jet : cs.inclusive_jets(ptmin)) {
                if (std::fabs(jet.eta()) < etamax) {
                    pass = true;
                    break;
                }
            }
        }
        if (pass) {
            for (const std::string& l : buf) out << l << "\n";
            nPassed++;
            sumwPassed += weight;
        }
        buf.clear();
        parts.clear();
        inEvent = false;
    };

    std::string line;
    while (std::getline(in, line)) {
        if (startsWith(line, "E ")) {
            flushEvent();
            inEvent = true;
            weight = eventWeight(line);
            buf.push_back(line);
            continue;
        }
        if (startsWith(line, "HepMC::")) {
            flushEvent();
            if (line.find("END_EVENT_LISTING") != std::string::npos) sawFooter = true;
            continue;
        }
        if (line.size() < 2 || !inEvent) continue;
        buf.push_back(line);
        if (startsWith(line, "P ")) {
            std::istringstream iss(line);
            std::string tag;
            long barcode;
            int pdgid, status;
            double px, py, pz, E, m;
            iss >> tag >> barcode >> pdgid >> px >> py >> pz >> E >> m >> status;
            if (status == 1 && !isExcludedLepton(pdgid)) {
                parts.emplace_back(px, py, pz, E);
            }
        }
    }

    if (!sawFooter) {
        std::cerr << "ERROR: " << inPath << " has no END_EVENT_LISTING footer "
                  << "(truncated JEWEL run?) -- refusing to filter it\n";
        out.close();
        std::remove(outPath.c_str());
        return 2;
    }

    out << "HepMC::IO_GenEvent-END_EVENT_LISTING\n\n";

    std::ofstream stats(outPath + ".stats");
    stats << "TOTAL " << nTotal << "\n"
          << "PASSED " << nPassed << "\n"
          << "SUMW_TOTAL " << sumwTotal << "\n"
          << "SUMW_PASSED " << sumwPassed << "\n";

    std::cout << "TOTAL " << nTotal << " PASSED " << nPassed
              << " SUMW_TOTAL " << sumwTotal << " SUMW_PASSED " << sumwPassed << "\n";
    return 0;
}
