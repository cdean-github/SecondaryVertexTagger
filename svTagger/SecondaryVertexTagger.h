// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SECONDARYVERTEXTAGGER_H
#define SECONDARYVERTEXTAGGER_H

#include "SecondaryVertexTagger_engine.h"
#include "SecondaryVertexTagger_nTuple.h"

#include <ffamodules/CDBInterface.h>
#include <fun4all/SubsysReco.h>
#include <phool/getClass.h>

//#include <KFParticle.h>

#include <TEntryList.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <filesystem>
#include <string>

class PHCompositeNode;
//class KFParticle;
// class TFile;

class SecondaryVertexTagger : public SubsysReco, public SecondaryVertexTagger_engine, public SecondaryVertexTagger_nTuple
{
public:
  SecondaryVertexTagger(const std::string &name = "SecondaryVertexTagger");

  ~SecondaryVertexTagger() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;
  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  //void Print(const std::string &what = "ALL") const override;

  // Use truth variables instead of reco
  void truthMode() { m_truth_mode = true; }
  void truthRecoMatch() { m_truth_reco_match = true; }

  void buildCharmDecays() { m_build_charm_decays = true; }

  void saveOutput() { m_save_output = true; }
  void setOutputName(const std::string &name) { m_outfile_name = name; }

  void magFieldFile(const std::string &name) { m_magField = name; }

  void setMaxTrackChi2nDoF(float value) { m_reco_chi2nDoF = value; }
  void setMaxTruthEta(float value) { m_max_eta = value; }
  void setMinTruthPt(float value) { m_min_pT = value; }

private:
  bool m_save_output = false;
  std::string m_outfile_name = "secondaryVertexInfo.root";
  TFile *m_outfile;

  void getField();

  std::string m_magField = "FIELDMAP_TRACKING";
};

#endif // SECONDARYVERTEXTAGGER_H
