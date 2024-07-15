#ifndef SECONDARYVERTEXTAGGER_ENGINE_H
#define SECONDARYVERTEXTAGGER_ENGINE_H

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <phool/getClass.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TDatabasePDG.h>

#include <KFVertex.h>

#include <KFParticle.h>
//#include <KFParticleBase.h>

#include <iostream>

//class KFParticle;
//class KFPVertex;
class PHCompositeNode;
class PHG4Particle;
class PHG4VtxPoint;
class SvtxTrackMap;
class SvtxTrack;

class SecondaryVertexTagger_engine
{
public:
protected:
  bool m_truth_mode = false;

  float m_reco_chi2nDoF = 10;
  float m_min_pT = 0.16;
  float m_max_eta = 1.1;

  void tagVertices(PHCompositeNode *topNode, std::vector<KFParticle>& tracks, std::vector<KFParticle>& vertices);

  void identify(const KFParticle &particle);

private:
  std::vector<KFParticle> makeAllParticles(PHCompositeNode *topNode);
  KFParticle makeParticlesFromReco();
  bool isGoodRecoParticle(KFParticle testParticle);
  KFParticle makeParticlesFromTruth();
  bool isGoodTruthParticle(KFParticle testParticle);

  std::vector<KFParticle> makeAllPrimaryVertices(PHCompositeNode *topNode);
  KFParticle makeVertexFromReco();
  KFParticle makeVertexFromTruth();

  float getParticleMass(const int PDGID);
  float getParticleCharge(const int PDGID);

  PHG4TruthInfoContainer *m_truthInfo{nullptr};
  PHG4Particle *m_g4_particle{nullptr};
  PHG4VtxPoint *m_g4_pv{nullptr};

  SvtxTrackMap *m_dst_trackmap{nullptr};
  SvtxTrack *m_dst_track{nullptr};

  SvtxVertexMap *m_dst_vertexmap{nullptr};
  SvtxVertex *m_dst_vertex{nullptr};
};

#endif
