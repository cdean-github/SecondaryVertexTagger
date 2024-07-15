/*****************/
/* Cameron Dean  */
/*   MIT 2024    */
/* cdean@bnl.gov */
/*****************/

#include "SecondaryVertexTagger_engine.h"
/*
 * Main body of tagger
 */
void SecondaryVertexTagger_engine::tagVertices(PHCompositeNode* topNode, std::vector<KFParticle>& tracks, std::vector<KFParticle>& vertices)
{
  tracks = makeAllParticles(topNode);
  vertices = makeAllPrimaryVertices(topNode);

  //Now lets see what KFParticle does for a PV reco
  //int nDaughters = tracks.size();
  //const KFParticle *trackArray[nDaughters];
  //for (unsigned int i = 0; i < tracks.size(); ++i)
  //{
  //  trackArray[i] = &tracks.at(i);
  //}
  //bool vtxFlag[nDaughters];
  //float chiCut = 3.5;

  //KFVertex kfPV;
  //kfPV.ConstructPrimaryVertex(trackArray, nDaughters, vtxFlag, chiCut);

  //std::cout << "The KFParticle reconstructed PV with a track chisq cut of " << chiCut << " has properties:" << std::endl;
  //identify(kfPV);
  
}

/*
 * Main loop to build particle object from reco or truth info
 */
std::vector<KFParticle> SecondaryVertexTagger_engine::makeAllParticles(PHCompositeNode *topNode)
{
  std::vector<KFParticle> daughterParticles;
  unsigned int trackID = 0;

  if (m_truth_mode)
  {
    m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetParticleRange();

    int trackableParticles[] = {11, 13, 211, 321, 2212};
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {
      m_g4_particle = iter->second;

      if (std::find(std::begin(trackableParticles), std::end(trackableParticles),
                    std::abs(m_g4_particle->get_pid())) != std::end(trackableParticles))
      {
        KFParticle particle = makeParticlesFromTruth();

        if (isGoodTruthParticle(particle))
        {
          daughterParticles.push_back(particle); /// Turn all dst tracks into KFP tracks
          daughterParticles[trackID].SetId(iter->first);
          ++trackID;
        }
      }
    }
  }
  else
  {
    m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    for (auto &iter : *m_dst_trackmap)
    {
      m_dst_track = iter.second;
      KFParticle particle = makeParticlesFromReco();

      if (isGoodRecoParticle(particle))
      {
        daughterParticles.push_back(particle); /// Turn all dst tracks into KFP tracks
        daughterParticles[trackID].SetId(iter.first);
        ++trackID;
      }
    }
  }

  return daughterParticles;
}

/*
 * Building particles from reconstructed info
 */
KFParticle SecondaryVertexTagger_engine::makeParticlesFromReco() /// Return a KFPTrack from track vector and covariance matrix. No mass or vertex constraints
{
  KFParticle kfp_particle;

  float f_trackParameters[6] = {m_dst_track->get_x(),
                                m_dst_track->get_y(),
                                m_dst_track->get_z(),
                                m_dst_track->get_px(),
                                m_dst_track->get_py(),
                                m_dst_track->get_pz()};

  float f_trackCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 6; ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
    {
      f_trackCovariance[iterate] = m_dst_track->get_error(i, j);
      ++iterate;
    }
  }

  kfp_particle.Create(f_trackParameters, f_trackCovariance, (Int_t) m_dst_track->get_charge(), -1);
  kfp_particle.NDF() = m_dst_track->get_ndf();
  kfp_particle.Chi2() = m_dst_track->get_chisq();
  kfp_particle.SetId(m_dst_track->get_id());

  return kfp_particle;
}

/*
 * Check to select reconstructed particles
 */
bool SecondaryVertexTagger_engine::isGoodRecoParticle(KFParticle testParticle) 
{
  float chi2nDoF = testParticle.GetChi2() / testParticle.GetNDF();
  return chi2nDoF <= m_reco_chi2nDoF ? true : false;
}

/*
 * Building particles from truth info
 */
KFParticle SecondaryVertexTagger_engine::makeParticlesFromTruth()
{
  KFParticle kfp_particle;
  PHG4VtxPoint *thisVtx = m_truthInfo->GetVtx(m_g4_particle->get_vtx_id());

  float f_trackParameters[6] = {(float) thisVtx->get_x(),
                                (float) thisVtx->get_y(),
                                (float) thisVtx->get_z(), 
                                (float) m_g4_particle->get_px(),
                                (float) m_g4_particle->get_py(),
                                (float) m_g4_particle->get_pz()};

  float f_trackCovariance[21] = {0.};

  int truthCharge = (Int_t) getParticleCharge(m_g4_particle->get_pid());
  float truthMass = getParticleMass(m_g4_particle->get_pid());

  kfp_particle.Create(f_trackParameters, f_trackCovariance, truthCharge, truthMass);
  kfp_particle.NDF() = 0;
  kfp_particle.Chi2() = 0;
  kfp_particle.SetId(m_g4_particle->get_track_id());
  kfp_particle.SetPDG(m_g4_particle->get_pid());

  return kfp_particle;
}

/*
 * Check to select truth particles
 */
bool SecondaryVertexTagger_engine::isGoodTruthParticle(KFParticle testParticle)
{
  bool isGood = false;
  float pT, pT_err, eta, eta_err;
  testParticle.GetPt(pT, pT_err);
  testParticle.GetEta(eta, eta_err);

  if ((pT >= m_min_pT) && (abs(eta) <= m_max_eta))
  {
    isGood = true;
  }

  return isGood;
}

/*
 * Main loop to build vertex object from reco or truth info
 */
std::vector<KFParticle> SecondaryVertexTagger_engine::makeAllPrimaryVertices(PHCompositeNode *topNode)
{
  unsigned int vertexID = 0;
  std::vector<KFParticle> primaryVertices;

  if (m_truth_mode)
  {
    m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::ConstVtxRange range = m_truthInfo->GetPrimaryVtxRange();

    for (PHG4TruthInfoContainer::ConstVtxIterator iter = range.first; iter != range.second; ++iter)
    {
      m_g4_pv = iter->second;
      primaryVertices.push_back(makeVertexFromTruth());
      primaryVertices[vertexID].SetId(m_g4_pv->get_id());
      ++vertexID;
    }
  }
  else
  {
    m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    auto globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!globalvertexmap)
    {
      std::cout << __FILE__ << ": no GlobalVertexMap" << std::endl;
      exit(1);
    }

    for (GlobalVertexMap::ConstIter iter = globalvertexmap->begin(); iter != globalvertexmap->end(); ++iter)
    {
      GlobalVertex *gvertex = iter->second;
      auto svtxiter = gvertex->find_vertexes(GlobalVertex::SVTX);
      // check that it contains a track vertex
      if (svtxiter == gvertex->end_vertexes())
      {
        continue;
      }

      auto svtxvertexvector = svtxiter->second;

      for (auto &vertex : svtxvertexvector)
      {
        m_dst_vertex = m_dst_vertexmap->find(vertex->get_id())->second;

        primaryVertices.push_back(makeVertexFromReco());
        primaryVertices[vertexID].SetId(gvertex->get_id());
        ++vertexID;
      }
    }
  }

  return primaryVertices;
}

/*
 * Building vertices from reconstructed info
 */
KFParticle SecondaryVertexTagger_engine::makeVertexFromReco()
{
  KFParticle kfp_vertex;

  float f_vertexParameters[6] = {m_dst_vertex->get_x(),
                                 m_dst_vertex->get_y(),
                                 m_dst_vertex->get_z(), 0, 0, 0};

  float f_vertexCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
    {
      f_vertexCovariance[iterate] = m_dst_vertex->get_error(i, j);
      ++iterate;
    }
  }

  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = m_dst_vertex->get_ndof();
  kfp_vertex.Chi2() = m_dst_vertex->get_chisq();

  return kfp_vertex;
}

/*
 * Building vertices from truth info
 */
KFParticle SecondaryVertexTagger_engine::makeVertexFromTruth()
{
  KFParticle kfp_vertex;

  float f_vertexParameters[6] = {(float) m_g4_pv->get_x(),
                                 (float) m_g4_pv->get_y(),
                                 (float) m_g4_pv->get_z(), 0, 0, 0};

  float f_vertexCovariance[21] = {0};

  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = 0;
  kfp_vertex.Chi2() = 0;

  return kfp_vertex;
}

/*
 * TParticlePDG call to get a particles mass
 */
float SecondaryVertexTagger_engine::getParticleMass(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->Mass();
}

/*
 * TParticlePDG call to get a particles charge
 */
float SecondaryVertexTagger_engine::getParticleCharge(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->Charge() / 3;
}

/*
 * General print for a KFParticle
 */
void SecondaryVertexTagger_engine::identify(const KFParticle &particle)
{
  std::cout << "Track ID: " << particle.Id() << std::endl;
  std::cout << "PDG ID: " << particle.GetPDG() << ", charge: " << (int)particle.GetQ() << ", mass: " << particle.GetMass() << " GeV" << std::endl;
  std::cout << "(px,py,pz) = (" << particle.GetPx() << " +/- " << std::sqrt(particle.GetCovariance(3, 3)) << ", ";
  std::cout <<                     particle.GetPy() << " +/- " << std::sqrt(particle.GetCovariance(4, 4)) << ", ";
  std::cout <<                     particle.GetPz() << " +/- " << std::sqrt(particle.GetCovariance(5, 5)) << ") GeV" << std::endl;
  std::cout << "(x,y,z) = (" << particle.GetX() << " +/- " << std::sqrt(particle.GetCovariance(0, 0)) << ", ";
  std::cout <<                  particle.GetY() << " +/- " << std::sqrt(particle.GetCovariance(1, 1)) << ", ";
  std::cout <<                  particle.GetZ() << " +/- " << std::sqrt(particle.GetCovariance(2, 2)) << ") cm\n" << std::endl;
}
