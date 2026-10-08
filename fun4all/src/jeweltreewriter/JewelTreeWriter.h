#ifndef JEWELTREEWRITER_H
#define JEWELTREEWRITER_H

// Flat-ntuple writer for the JEWEL calo-image pipeline.
//
// Reads calibrated calorimeter towers and clustered jets off the node tree and writes three
// flat TTrees into its own output file:
//   towertree  -- per-event EMCal(retowered)/IHCal/OHCal tower energies on the fixed
//                 eta x phi grid, plus the (constant) etabin/phibin index arrays
//   jettree    -- reco and truth jet pt/eta/phi above a pt threshold, each with the
//                 kinematics of its individual constituents. Truth constituents additionally
//                 carry energy and PDG id (tjet_constit_e / tjet_constit_pid), which is what
//                 makes their rapidity recoverable: FastJet clusters in y, but eta is what is
//                 stored, and the two differ for a massive particle at low pt. Reco
//                 constituents carry the Jet::SRC they came from (rjet_constit_src), i.e.
//                 which calorimeter read the tower out, so the jet can be split by layer.
//   globaltree -- the event vertex
//
// The branch names are fixed by what the downstream calo-image and QA macros already read
// (emcalenergy, rjet_pt, vx, ...) -- renaming them breaks draw_event_v2.C and QA_JewelTrees.C.
//
// Constituent kinematics need explaining: JetContainer stores each jet's constituents only as
// (Jet::SRC, index) pairs pointing back into the *original* tower/particle containers, with no
// kinematics attached. To turn those back into pt/eta/phi this module re-runs the same
// JetInput objects the jets were clustered from (see add_tower_constituent_source() /
// set_truth_constituent_embed_flag(), whose defaults mirror the standard sPHENIX tower and
// truth jet reco), which hands back single-particle Jets carrying both the (SRC, index) key
// and the kinematics. Those are matched against each jet's comp_vec entries. End() reports
// what fraction resolved; anything below ~100% means the configured inputs here don't match
// the ones the jets were actually clustered from.

#include <fun4all/SubsysReco.h>

#include <jetbase/Jet.h>

#include <map>
#include <string>
#include <vector>

class JetInput;
class PHCompositeNode;
class TFile;
class TTree;

class JewelTreeWriter : public SubsysReco
{
 public:
  explicit JewelTreeWriter(const std::string &name = "JewelTreeWriter",
                           const std::string &outfile = "jewel_trees.root");

  ~JewelTreeWriter() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  //! output file, if not given in the constructor
  void set_output_file(const std::string &outfile) { m_outfilename = outfile; }

  //! jets below this pt are dropped from jettree entirely (not recoverable downstream)
  void set_jet_pt_min(float ptmin) { m_jet_pt_min = ptmin; }

  //! calorimeter tower nodes written to towertree
  void set_emcal_tower_node(const std::string &node) { m_emcal_node = node; }
  void set_ihcal_tower_node(const std::string &node) { m_ihcal_node = node; }
  void set_ohcal_tower_node(const std::string &node) { m_ohcal_node = node; }

  //! tower grid the three energy arrays are flattened onto, as ieta * nphi + iphi
  void set_tower_grid(int neta, int nphi)
  {
    m_neta = neta;
    m_nphi = nphi;
  }

  //! jet containers written to jettree
  void set_reco_jet_node(const std::string &node) { m_reco_jet_node = node; }
  void set_truth_jet_node(const std::string &node) { m_truth_jet_node = node; }

  //! vertex node written to globaltree
  void set_vertex_node(const std::string &node) { m_vertex_node = node; }

  //! Replace the default tower constituent sources (IHCal + OHCal + retowered EMCal). Call
  //! once per source; the first call clears the defaults. These must match the inputs the
  //! reco jets were clustered from, or their constituents won't resolve.
  void add_tower_constituent_source(Jet::SRC src);

  //! node prefix handed to the tower JetInputs (default "TOWERINFO_CALIB")
  void set_tower_node_prefix(const std::string &prefix) { m_tower_prefix = prefix; }

  //! embedding id of the truth particles the truth jets were clustered from
  void set_truth_constituent_embed_flag(int flag) { m_truth_embed_flag = flag; }

  //! skip the constituent branches entirely (the *_constit_* branches are then empty)
  void set_do_constituents(bool doit) { m_do_constituents = doit; }

 private:
  void FillTowers(PHCompositeNode *topNode, const std::string &nodename, std::vector<float> &out);

  //! fills jet kinematics (pt >= m_jet_pt_min only) plus, per jet, the kinematics of every
  //! constituent resolvable through m_constit_map (rebuilt each event in process_event).
  //! ce and cpid are optional and only meaningful for truth jets; pass nullptr for towers.
  //! csrc is optional and only meaningful for tower jets, where the Jet::SRC of a constituent
  //! names the calorimeter it was read out of; truth constituents are all Jet::PARTICLE.
  void FillJets(PHCompositeNode *topNode, const std::string &nodename,
                std::vector<float> &pt, std::vector<float> &eta, std::vector<float> &phi,
                std::vector<std::vector<float>> &cpt,
                std::vector<std::vector<float>> &ceta,
                std::vector<std::vector<float>> &cphi,
                std::vector<std::vector<float>> *ce = nullptr,
                std::vector<std::vector<int>> *cpid = nullptr,
                std::vector<std::vector<int>> *csrc = nullptr);

  //! re-runs the configured JetInputs and indexes the single-particle Jets they return by
  //! their own (Jet::SRC, index) key, so jet comp_vec entries can be looked up
  void BuildConstituentMap(PHCompositeNode *topNode);
  void ClearConstituentMap();

  std::string m_outfilename;

  std::string m_emcal_node{"TOWERINFO_CALIB_CEMC_RETOWER"};
  std::string m_ihcal_node{"TOWERINFO_CALIB_HCALIN"};
  std::string m_ohcal_node{"TOWERINFO_CALIB_HCALOUT"};
  std::string m_reco_jet_node{"AntiKt_Tower_r04"};
  std::string m_truth_jet_node{"AntiKt_TruthFixed_r04"};
  std::string m_vertex_node{"GlobalVertexMap"};
  std::string m_tower_prefix{"TOWERINFO_CALIB"};

  int m_neta{24};
  int m_nphi{64};
  float m_jet_pt_min{10.};
  int m_truth_embed_flag{0};
  bool m_do_constituents{true};

  std::vector<Jet::SRC> m_tower_constituent_srcs{Jet::HCALIN_TOWERINFO,
                                                 Jet::HCALOUT_TOWERINFO,
                                                 Jet::CEMC_TOWERINFO_RETOWER};
  bool m_tower_srcs_overridden{false};

  //! owned; built in Init(), deleted in the destructor
  std::vector<JetInput *> m_constituent_inputs;
  //! owned per event: the single-particle Jets the inputs hand back, cleared after each event
  std::vector<Jet *> m_constituent_jets;
  //! non-owning index into m_constituent_jets, keyed by each seed's own (SRC, index)
  std::map<Jet::TYPE_comp, Jet *> m_constit_map;

  TFile *m_outfile{nullptr};
  TTree *m_towertree{nullptr};
  TTree *m_jettree{nullptr};
  TTree *m_globaltree{nullptr};

  std::vector<float> m_emcalenergy;
  std::vector<float> m_ihcalenergy;
  std::vector<float> m_ohcalenergy;
  std::vector<int> m_etabin;
  std::vector<int> m_phibin;

  std::vector<float> m_rjet_pt;
  std::vector<float> m_rjet_eta;
  std::vector<float> m_rjet_phi;
  std::vector<std::vector<float>> m_rjet_constit_pt;
  std::vector<std::vector<float>> m_rjet_constit_eta;
  std::vector<std::vector<float>> m_rjet_constit_phi;
  //! Jet::SRC of each reco constituent, cast to int: which calorimeter read that tower out
  std::vector<std::vector<int>> m_rjet_constit_src;

  std::vector<float> m_tjet_pt;
  std::vector<float> m_tjet_eta;
  std::vector<float> m_tjet_phi;
  std::vector<std::vector<float>> m_tjet_constit_pt;
  std::vector<std::vector<float>> m_tjet_constit_eta;
  std::vector<std::vector<float>> m_tjet_constit_phi;
  std::vector<std::vector<float>> m_tjet_constit_e;
  std::vector<std::vector<int>> m_tjet_constit_pid;

  float m_vx{-999.};
  float m_vy{-999.};
  float m_vz{-999.};

  long long m_ncomp_total{0};
  long long m_ncomp_matched{0};
  long long m_npid_total{0};
  long long m_npid_matched{0};
};

#endif  // JEWELTREEWRITER_H
