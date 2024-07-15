#ifndef SECONDARYVERTEXTAGGER_NTUPLE_H
#define SECONDARYVERTEXTAGGER_NTUPLE_H

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase_historic/SvtxPHG4ParticleMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <KFParticle.h>

#include <TLorentzVector.h>

#include <limits>
#include <string>  // for string
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float>>+;
#endif

class PHCompositeNode;
class PHG4Particle;
class TTree;

class SecondaryVertexTagger_nTuple
{
 public:
  /// Constructor
  SecondaryVertexTagger_nTuple();

  /// Destructor
  ~SecondaryVertexTagger_nTuple();

  /// Initialises required branches based off the user selection (number of tracks, PV constraints etc ) and sets branch names if specified
  void initializeBranches();

  void fillBranch(PHCompositeNode* topNode, std::vector<KFParticle> tracks, std::vector<KFParticle> vertices, bool isTruth);

  void clearVectors();

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  PHG4Particle *getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode);

  void buildD0Decays();

 protected:

  bool m_truth_reco_match = false;
  bool m_build_charm_decays = false;

 private:
  TTree *m_tree = nullptr;

  std::vector<float> PV_x;
  std::vector<float> PV_y;
  std::vector<float> PV_z;
  std::vector<float> particle_x;
  std::vector<float> particle_y;
  std::vector<float> particle_z;
  std::vector<int> particle_PID;
  std::vector<int> particle_charge;
  std::vector<float> particle_px;
  std::vector<float> particle_py;
  std::vector<float> particle_pz;
  std::vector<float> particle_pT;
  std::vector<float> particle_one_over_pT;
  std::vector<float> particle_eta;
  std::vector<float> particle_phi;
  std::vector<std::vector<float>> particle_hit_x;
  std::vector<std::vector<float>> particle_hit_y;
  std::vector<std::vector<float>> particle_hit_z;

  std::vector<float> particle_true_px;
  std::vector<float> particle_true_py;
  std::vector<float> particle_true_pz;
  std::vector<float> particle_true_pT;
  std::vector<float> particle_true_one_over_pT;

  std::vector<float> m_KminusPiplus;
  std::vector<float> m_KplusPiminus;

  PHG4TruthInfoContainer *m_truthInfo{nullptr};
  SvtxEvalStack *m_svtx_evalstack{nullptr};
  SvtxTruthEval *trutheval{nullptr};
  SvtxTrackEval *trackeval{nullptr};
  SvtxTrackMap *m_dst_trackmap{nullptr};
  SvtxTrack *m_dst_track{nullptr};
  TrkrClusterContainer *trktClusterContainer{nullptr};
  ActsGeometry *actsGeom{nullptr};
};

#endif
