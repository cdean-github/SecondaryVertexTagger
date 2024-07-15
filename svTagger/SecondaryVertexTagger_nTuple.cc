#include "SecondaryVertexTagger_nTuple.h"

#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <KFParticle.h>

#include <Rtypes.h>
#include <TString.h>  // for TString, operator+
#include <TTree.h>

#include <algorithm>  // for max
#include <cmath>
#include <cstdlib>  // for abs, size_t
#include <map>      // for map, _Rb_tree_iterator, map<>:...

class PHCompositeNode;
class PHNode;

SecondaryVertexTagger_nTuple::SecondaryVertexTagger_nTuple() {} // Constructor

SecondaryVertexTagger_nTuple::~SecondaryVertexTagger_nTuple() {}

void SecondaryVertexTagger_nTuple::initializeBranches()
{
  delete m_tree;
  m_tree = new TTree("DecayTree", "DecayTree");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);  // Save the output file every 5MB

  m_tree->Branch("primary_vertex_x", &PV_x);
  m_tree->Branch("primary_vertex_y", &PV_y);
  m_tree->Branch("primary_vertex_z", &PV_z);

  m_tree->Branch("particle_production_x", &particle_x);
  m_tree->Branch("particle_production_y", &particle_y);
  m_tree->Branch("particle_production_z", &particle_z);
  m_tree->Branch("particle_PID", &particle_PID);
  m_tree->Branch("particle_charge", &particle_charge);
  m_tree->Branch("particle_px", &particle_px);
  m_tree->Branch("particle_py", &particle_py);
  m_tree->Branch("particle_pz", &particle_pz);
  m_tree->Branch("particle_pT", &particle_pT);
  m_tree->Branch("particle_one_over_pT", &particle_one_over_pT);
  m_tree->Branch("particle_eta", &particle_eta);
  m_tree->Branch("particle_phi", &particle_phi);
  m_tree->Branch("particle_hit_x", &particle_hit_x);
  m_tree->Branch("particle_hit_y", &particle_hit_y);
  m_tree->Branch("particle_hit_z", &particle_hit_z);

  if (m_truth_reco_match)
  {
    m_tree->Branch("particle_true_px", &particle_true_px);
    m_tree->Branch("particle_true_py", &particle_true_py);
    m_tree->Branch("particle_true_pz", &particle_true_pz);
    m_tree->Branch("particle_true_pT", &particle_true_pT);
    m_tree->Branch("particle_true_one_over_pT", &particle_true_one_over_pT);
  }

  if (m_build_charm_decays)
  {
    m_tree->Branch("mass_KminusPiplus", &m_KminusPiplus);
    m_tree->Branch("mass_KplusPiminus", &m_KplusPiminus);
  }
}

void SecondaryVertexTagger_nTuple::fillBranch(PHCompositeNode* topNode, std::vector<KFParticle> tracks, std::vector<KFParticle> vertices, bool isTruth) 
{ 
  if (!isTruth)
  {
    PHNodeIterator nodeIter(topNode);
    PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxTrackMap"));
    if (findNode)
    {
      m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    }

    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("TRKR_CLUSTER"));
    if (findNode)
    {
      trktClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }

    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("ActsGeometry"));
    if (findNode)
    {
      actsGeom = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    }
  }

  for (auto vertex : vertices)
  {
    PV_x.push_back(vertex.GetX());
    PV_y.push_back(vertex.GetY());
    PV_z.push_back(vertex.GetZ());
  }
 
  for (auto track : tracks)
  {
    int pdgID = std::numeric_limits<int>::quiet_NaN();
    PHG4Particle* truthMatchedParticle{nullptr};
    if (m_truth_reco_match)
    {
      m_dst_track = getTrack(track.Id(), m_dst_trackmap);
      truthMatchedParticle = getTruthTrack(m_dst_track, topNode);
      bool isParticleValid = truthMatchedParticle == nullptr ? false : true;

      pdgID = isParticleValid ? truthMatchedParticle->get_pid() : std::numeric_limits<int>::quiet_NaN();
      float true_px = isParticleValid ? truthMatchedParticle->get_px() : std::numeric_limits<float>::quiet_NaN();
      float true_py = isParticleValid ? truthMatchedParticle->get_py() : std::numeric_limits<float>::quiet_NaN();
      float true_pz = isParticleValid ? truthMatchedParticle->get_pz() : std::numeric_limits<float>::quiet_NaN();
      float true_pT = isParticleValid ? sqrt(pow(true_px, 2) + pow(true_py, 2)) : std::numeric_limits<float>::quiet_NaN();
      float true_one_over_pT = isParticleValid ? 1./true_pT : std::numeric_limits<float>::quiet_NaN();

      particle_true_px.push_back(true_px);
      particle_true_py.push_back(true_py);
      particle_true_pz.push_back(true_pz);
      particle_true_pT.push_back(true_pT);
      particle_true_one_over_pT.push_back(true_one_over_pT);
    }
    
    particle_x.push_back(track.GetX());
    particle_y.push_back(track.GetY());
    particle_z.push_back(track.GetZ());
    particle_PID.push_back(pdgID);
    particle_charge.push_back((int) track.GetQ());
    particle_px.push_back(track.GetPx());
    particle_py.push_back(track.GetPy());
    particle_pz.push_back(track.GetPz());

    float pT, pT_err, eta, eta_err, phi, phi_err;
    track.GetPt(pT, pT_err);
    track.GetEta(eta, eta_err);
    track.GetPhi(phi, phi_err);

    particle_pT.push_back(pT);
    particle_one_over_pT.push_back(1/pT);
    particle_eta.push_back(eta);
    particle_phi.push_back(phi);

    std::vector<float> x, y, z;

    if (isTruth)
    {
      if (!m_svtx_evalstack)
      {
        m_svtx_evalstack = new SvtxEvalStack(topNode);
        trutheval = m_svtx_evalstack->get_truth_eval();
      }

      m_svtx_evalstack->next_event(topNode);

      PHG4Particle* truthParticle = trutheval->get_particle(track.Id());
      std::set<PHG4Hit*> truthHits = trutheval->all_truth_hits(truthParticle);

      for (auto hit : truthHits)
      {
        if (hit->get_layer() <= 6)
        {
          x.push_back(hit->get_x(0));
          y.push_back(hit->get_y(0));
          z.push_back(hit->get_z(0)); 
        }
      }
      if (x.size() == 0)
      {
        x.push_back(std::numeric_limits<float>::quiet_NaN());
        x.push_back(std::numeric_limits<float>::quiet_NaN());
        z.push_back(std::numeric_limits<float>::quiet_NaN());
      }
    }
    else
    {
      m_dst_track = getTrack(track.Id(), m_dst_trackmap);
      TrackSeed *silseed = m_dst_track->get_silicon_seed();

      if (silseed)
      {
        for (auto cluster_iter = silseed->begin_cluster_keys(); cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
        {
          const auto &stateckey = *cluster_iter;
          TrkrCluster *cluster = trktClusterContainer->findCluster(stateckey);
          Acts::Vector3 global = actsGeom->getGlobalPosition(stateckey, cluster);
          x.push_back(global.x());
          y.push_back(global.y());
          z.push_back(global.z());
        }
      }
    }

    particle_hit_x.push_back(x);
    particle_hit_y.push_back(y);
    particle_hit_z.push_back(z);
  } 

  if (m_build_charm_decays) buildD0Decays(); 

  m_tree->Fill();
  clearVectors();
}

void SecondaryVertexTagger_nTuple::clearVectors()
{
  PV_x.clear();
  PV_y.clear();
  PV_z.clear();
  particle_x.clear();
  particle_y.clear();
  particle_z.clear();
  particle_PID.clear();
  particle_charge.clear();
  particle_px.clear();
  particle_py.clear();
  particle_pz.clear();
  particle_pT.clear();
  particle_one_over_pT.clear();
  particle_eta.clear();
  particle_phi.clear();
  particle_hit_x.clear();
  particle_hit_y.clear();
  particle_hit_z.clear();
  particle_true_px.clear();
  particle_true_py.clear();
  particle_true_pz.clear();
  particle_true_pT.clear();
  particle_true_one_over_pT.clear();
  m_KminusPiplus.clear();
  m_KplusPiminus.clear();
}

SvtxTrack *SecondaryVertexTagger_nTuple::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = nullptr;

  for (auto &iter : *trackmap)
  {
    if (iter.first == track_id)
    {
      matched_track = iter.second;
    }
  }

  return matched_track;
}

PHG4Particle *SecondaryVertexTagger_nTuple::getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode)
{
  /*
   * There are two methods for getting the truth rack from the reco track
   * 1. (recommended) Use the reco -> truth tables (requires SvtxPHG4ParticleMap). Introduced Summer of 2022
   * 2. Get truth track via nClusters. Older method and will work with older DSTs
   */

  PHG4Particle *particle = nullptr;

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxPHG4ParticleMap"));
  if (findNode)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "SecondaryVertexTagger truth matching: G4TruthInfo does not exist" << std::endl;
    }

    SvtxPHG4ParticleMap *dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");

    std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(thisTrack->get_id());
    const auto &best_weight = truth_set.rbegin();
    int best_truth_id = *best_weight->second.rbegin();
    particle = m_truthInfo->GetParticle(best_truth_id);
  }
  else
  {
    std::cout << __FILE__ << ": SvtxPHG4ParticleMap not found, reverting to max_truth_particle_by_nclusters()" << std::endl;

    if (!m_svtx_evalstack)
    {
      m_svtx_evalstack = new SvtxEvalStack(topNode);
      trackeval = m_svtx_evalstack->get_track_eval();
    }

    m_svtx_evalstack->next_event(topNode);

    particle = trackeval->max_truth_particle_by_nclusters(thisTrack);
  }
  return particle;
}

void SecondaryVertexTagger_nTuple::buildD0Decays()
{
  float pionMass = 0.13957; //GeV
  float kaonMass = 0.49368; //GeV

  unsigned int nTracks = particle_charge.size();
  for (unsigned int i = 0; i < nTracks - 1; ++i)
  {
    for (unsigned int j = i + 1; j < nTracks; ++j)
    {
      if (particle_charge.at(i) == particle_charge.at(j)) continue;

      float track_1_px = particle_px.at(i);
      float track_1_py = particle_py.at(i);
      float track_1_pz = particle_pz.at(i);

      float track_2_px = particle_px.at(j);
      float track_2_py = particle_py.at(j);
      float track_2_pz = particle_pz.at(j);

      float track_1_energy_pion = sqrt(pow(track_1_px, 2) + pow(track_1_py, 2) + pow(track_1_pz, 2) + pow(pionMass, 2));
      float track_1_energy_kaon = sqrt(pow(track_1_px, 2) + pow(track_1_py, 2) + pow(track_1_pz, 2) + pow(kaonMass, 2));
      float track_2_energy_pion = sqrt(pow(track_2_px, 2) + pow(track_2_py, 2) + pow(track_2_pz, 2) + pow(pionMass, 2));
      float track_2_energy_kaon = sqrt(pow(track_2_px, 2) + pow(track_2_py, 2) + pow(track_2_pz, 2) + pow(kaonMass, 2));

      TLorentzVector vector1_pion(track_1_px, track_1_py, track_1_pz, track_1_energy_pion);
      TLorentzVector vector1_kaon(track_1_px, track_1_py, track_1_pz, track_1_energy_kaon);
      TLorentzVector vector2_pion(track_2_px, track_2_py, track_2_pz, track_2_energy_pion);
      TLorentzVector vector2_kaon(track_2_px, track_2_py, track_2_pz, track_2_energy_kaon);

      TLorentzVector vectorSum_Kminus_piplus = particle_charge.at(i) == -1 ? vector1_kaon + vector2_pion : vector1_pion + vector2_kaon;
      TLorentzVector vectorSum_Kplus_piminus = particle_charge.at(i) == +1 ? vector1_kaon + vector2_pion : vector1_pion + vector2_kaon;

      float mass_Kminus_piplus = vectorSum_Kminus_piplus.M();
      float mass_Kplus_piminus = vectorSum_Kplus_piminus.M();

      float massRange[2] = {1.7, 2.1};
      if ((massRange[0] <= mass_Kminus_piplus) && (mass_Kminus_piplus <= massRange[1])) m_KminusPiplus.push_back(mass_Kminus_piplus); 
      if ((massRange[0] <= mass_Kplus_piminus) && (mass_Kplus_piminus <= massRange[1])) m_KplusPiminus.push_back(mass_Kplus_piminus); 
    }
  } 
}
