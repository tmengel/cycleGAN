#include "JewelTreeWriter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetInput.h>
#include <jetbase/TowerJetInput.h>

#include <g4jets/TruthJetInput.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <utility>

JewelTreeWriter::JewelTreeWriter(const std::string &name, const std::string &outfile)
  : SubsysReco(name)
  , m_outfilename(outfile)
{
}

JewelTreeWriter::~JewelTreeWriter()
{
  ClearConstituentMap();
  for (JetInput *input : m_constituent_inputs)
  {
    delete input;
  }
}

void JewelTreeWriter::add_tower_constituent_source(Jet::SRC src)
{
  if (!m_tower_srcs_overridden)
  {
    m_tower_constituent_srcs.clear();
    m_tower_srcs_overridden = true;
  }
  m_tower_constituent_srcs.push_back(src);
}

int JewelTreeWriter::Init(PHCompositeNode * /*topNode*/)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");

  m_towertree = new TTree("towertree", "JEWEL calo towers");
  m_towertree->Branch("emcalenergy", &m_emcalenergy);
  m_towertree->Branch("ihcalenergy", &m_ihcalenergy);
  m_towertree->Branch("ohcalenergy", &m_ohcalenergy);
  m_towertree->Branch("phibin", &m_phibin);
  m_towertree->Branch("etabin", &m_etabin);

  m_jettree = new TTree("jettree", "JEWEL truth+reco jets");
  m_jettree->Branch("rjet_pt", &m_rjet_pt);
  m_jettree->Branch("rjet_eta", &m_rjet_eta);
  m_jettree->Branch("rjet_phi", &m_rjet_phi);
  m_jettree->Branch("rjet_constit_pt", &m_rjet_constit_pt);
  m_jettree->Branch("rjet_constit_eta", &m_rjet_constit_eta);
  m_jettree->Branch("rjet_constit_phi", &m_rjet_constit_phi);
  m_jettree->Branch("rjet_constit_src", &m_rjet_constit_src);
  m_jettree->Branch("tjet_pt", &m_tjet_pt);
  m_jettree->Branch("tjet_eta", &m_tjet_eta);
  m_jettree->Branch("tjet_phi", &m_tjet_phi);
  m_jettree->Branch("tjet_constit_pt", &m_tjet_constit_pt);
  m_jettree->Branch("tjet_constit_eta", &m_tjet_constit_eta);
  m_jettree->Branch("tjet_constit_phi", &m_tjet_constit_phi);
  m_jettree->Branch("tjet_constit_e", &m_tjet_constit_e);
  m_jettree->Branch("tjet_constit_pid", &m_tjet_constit_pid);

  m_globaltree = new TTree("globaltree", "JEWEL event vertex");
  m_globaltree->Branch("vx", &m_vx);
  m_globaltree->Branch("vy", &m_vy);
  m_globaltree->Branch("vz", &m_vz);

  if (m_do_constituents)
  {
    for (const Jet::SRC src : m_tower_constituent_srcs)
    {
      m_constituent_inputs.push_back(new TowerJetInput(src, m_tower_prefix));
    }
    TruthJetInput *truthinput = new TruthJetInput(Jet::PARTICLE);
    truthinput->add_embedding_flag(m_truth_embed_flag);
    m_constituent_inputs.push_back(truthinput);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JewelTreeWriter::process_event(PHCompositeNode *topNode)
{
  FillTowers(topNode, m_emcal_node, m_emcalenergy);
  FillTowers(topNode, m_ihcal_node, m_ihcalenergy);
  FillTowers(topNode, m_ohcal_node, m_ohcalenergy);

  // the bin index arrays are the same grid every event -- fill once, on the first event
  if (m_phibin.empty())
  {
    for (int ieta = 0; ieta < m_neta; ieta++)
    {
      for (int iphi = 0; iphi < m_nphi; iphi++)
      {
        m_etabin.push_back(ieta);
        m_phibin.push_back(iphi);
      }
    }
  }

  BuildConstituentMap(topNode);

  FillJets(topNode, m_truth_jet_node,
           m_tjet_pt, m_tjet_eta, m_tjet_phi,
           m_tjet_constit_pt, m_tjet_constit_eta, m_tjet_constit_phi,
           &m_tjet_constit_e, &m_tjet_constit_pid);
  FillJets(topNode, m_reco_jet_node,
           m_rjet_pt, m_rjet_eta, m_rjet_phi,
           m_rjet_constit_pt, m_rjet_constit_eta, m_rjet_constit_phi,
           nullptr, nullptr, &m_rjet_constit_src);

  ClearConstituentMap();

  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, m_vertex_node);
  GlobalVertex *vtx = (vertexmap && !vertexmap->empty()) ? vertexmap->begin()->second : nullptr;
  m_vx = vtx ? vtx->get_x() : -999.;
  m_vy = vtx ? vtx->get_y() : -999.;
  m_vz = vtx ? vtx->get_z() : -999.;

  m_towertree->Fill();
  m_jettree->Fill();
  m_globaltree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int JewelTreeWriter::End(PHCompositeNode * /*topNode*/)
{
  m_outfile->cd();
  m_towertree->Write();
  m_jettree->Write();
  m_globaltree->Write();
  m_outfile->Close();

  if (m_ncomp_total > 0)
  {
    std::cout << Name() << ": resolved " << m_ncomp_matched << "/" << m_ncomp_total
              << " jet constituent entries ("
              << (100. * static_cast<double>(m_ncomp_matched) / static_cast<double>(m_ncomp_total))
              << "%) via the rebuilt seed-particle lookup." << std::endl;
  }

  if (m_npid_total > 0)
  {
    std::cout << Name() << ": resolved " << m_npid_matched << "/" << m_npid_total
              << " truth constituent PDG ids ("
              << (100. * static_cast<double>(m_npid_matched) / static_cast<double>(m_npid_total))
              << "%). Well below 100% means the truth comp key is not the G4 track id."
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void JewelTreeWriter::FillTowers(PHCompositeNode *topNode, const std::string &nodename,
                                 std::vector<float> &out)
{
  out.clear();
  out.resize(static_cast<size_t>(m_neta * m_nphi), 0.);

  TowerInfoContainer *towers = findNode::getClass<TowerInfoContainer>(topNode, nodename);
  if (!towers)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing tower node " << nodename << std::endl;
    }
    return;
  }

  const size_t ntowers = towers->size();
  for (size_t channel = 0; channel < ntowers; channel++)
  {
    const unsigned int towerkey = towers->encode_key(static_cast<unsigned int>(channel));
    const int ieta = static_cast<int>(towers->getTowerEtaBin(towerkey));
    const int iphi = static_cast<int>(towers->getTowerPhiBin(towerkey));
    if (ieta < 0 || ieta >= m_neta || iphi < 0 || iphi >= m_nphi)
    {
      continue;
    }
    out[static_cast<size_t>(ieta * m_nphi + iphi)] =
        towers->get_tower_at_channel(static_cast<int>(channel))->get_energy();
  }
}

void JewelTreeWriter::BuildConstituentMap(PHCompositeNode *topNode)
{
  if (!m_do_constituents)
  {
    return;
  }

  for (JetInput *input : m_constituent_inputs)
  {
    std::vector<Jet *> seeds = input->get_input(topNode);
    for (Jet *seed : seeds)
    {
      // each seed carries exactly one (SRC, index) entry -- the key the clustered jets
      // reference their constituents by
      if (!seed->get_comp_vec().empty())
      {
        m_constit_map[seed->get_comp_vec().front()] = seed;
      }
      m_constituent_jets.push_back(seed);
    }
  }
}

void JewelTreeWriter::ClearConstituentMap()
{
  m_constit_map.clear();
  for (Jet *seed : m_constituent_jets)
  {
    delete seed;
  }
  m_constituent_jets.clear();
}

void JewelTreeWriter::FillJets(PHCompositeNode *topNode, const std::string &nodename,
                               std::vector<float> &pt, std::vector<float> &eta,
                               std::vector<float> &phi,
                               std::vector<std::vector<float>> &cpt,
                               std::vector<std::vector<float>> &ceta,
                               std::vector<std::vector<float>> &cphi,
                               std::vector<std::vector<float>> *ce,
                               std::vector<std::vector<int>> *cpid,
                               std::vector<std::vector<int>> *csrc)
{
  pt.clear();
  eta.clear();
  phi.clear();
  cpt.clear();
  ceta.clear();
  cphi.clear();
  if (ce) ce->clear();
  if (cpid) cpid->clear();
  if (csrc) csrc->clear();

  // Only needed for the PDG ids. The (SRC, index) key a truth seed carries is assumed to be the
  // G4 track id; rather than trust that, every lookup is checked against the seed's momentum and
  // End() reports how many resolved, the same way the constituent match rate is reported.
  PHG4TruthInfoContainer *truthinfo =
      cpid ? findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo") : nullptr;

  JetContainer *jets = findNode::getClass<JetContainer>(topNode, nodename);
  if (!jets)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing jet node " << nodename << std::endl;
    }
    return;
  }

  for (auto *jet : *jets)
  {
    if (jet->get_pt() < m_jet_pt_min)
    {
      continue;
    }
    pt.push_back(jet->get_pt());
    eta.push_back(jet->get_eta());
    phi.push_back(jet->get_phi());

    std::vector<float> jet_cpt;
    std::vector<float> jet_ceta;
    std::vector<float> jet_cphi;
    std::vector<float> jet_ce;
    std::vector<int> jet_cpid;
    std::vector<int> jet_csrc;
    if (m_do_constituents)
    {
      for (const auto &comp : jet->get_comp_vec())
      {
        m_ncomp_total++;
        auto it = m_constit_map.find(comp);
        if (it == m_constit_map.end())
        {
          continue;
        }
        m_ncomp_matched++;
        jet_cpt.push_back(it->second->get_pt());
        jet_ceta.push_back(it->second->get_eta());
        jet_cphi.push_back(it->second->get_phi());
        if (csrc)
        {
          // comp.first is the container the jet took this constituent from, i.e. the
          // calorimeter that read the tower out. Kept as the raw Jet::SRC value so it stays
          // meaningful if the tower sources are reconfigured; see add_tower_constituent_source().
          jet_csrc.push_back(static_cast<int>(comp.first));
        }
        if (ce)
        {
          jet_ce.push_back(it->second->get_e());
        }
        if (cpid)
        {
          int pid = 0;
          PHG4Particle *part =
              truthinfo ? truthinfo->GetParticle(static_cast<int>(comp.second)) : nullptr;
          if (part)
          {
            const double dpx = part->get_px() - it->second->get_px();
            const double dpy = part->get_py() - it->second->get_py();
            const double dpz = part->get_pz() - it->second->get_pz();
            // same particle, same numbers: anything else means the key is not the track id
            if (std::sqrt(dpx * dpx + dpy * dpy + dpz * dpz) < 1e-3)
            {
              pid = part->get_pid();
            }
          }
          m_npid_total++;
          if (pid != 0)
          {
            m_npid_matched++;
          }
          jet_cpid.push_back(pid);
        }
      }
    }
    cpt.push_back(jet_cpt);
    ceta.push_back(jet_ceta);
    cphi.push_back(jet_cphi);
    if (ce) ce->push_back(jet_ce);
    if (cpid) cpid->push_back(jet_cpid);
    if (csrc) csrc->push_back(jet_csrc);
  }
}
