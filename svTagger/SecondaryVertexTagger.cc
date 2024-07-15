//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in SecondaryVertexTagger.h.
//
// SecondaryVertexTagger(const std::string &name = "SecondaryVertexTagger")
// everything is keyed to SecondaryVertexTagger, duplicate names do work but it
// makes e.g. finding culprits in logs difficult or getting a pointer to the
// module from the command line
//
// SecondaryVertexTagger::~SecondaryVertexTagger()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the
// node tree
//
// int SecondaryVertexTagger::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer.
// You can create historgrams here or put objects on the node tree but be aware
// that modules which haven't been registered yet did not put antyhing on the
// node tree
//
// int SecondaryVertexTagger::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's
// action depends on what else is around. Last chance to put nodes under the DST
// Node We mix events during readback if branches are added after the first
// event
//
// int SecondaryVertexTagger::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT;
//   proceed but do not save this event in output (needs output manager
//   setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT;
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int SecondaryVertexTagger::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs
// clearing after each event, this is the place to do that. The nodes under the
// DST node are cleared by the framework
//
// int SecondaryVertexTagger::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int SecondaryVertexTagger::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int SecondaryVertexTagger::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void SecondaryVertexTagger::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

/*****************/
/* Cameron Dean  */
/*   MIT 2024    */
/* cdean@bnl.gov */
/*****************/

#include "SecondaryVertexTagger.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

int candidateCounter = 0;

//____________________________________________________________________________..
SecondaryVertexTagger::SecondaryVertexTagger(const std::string &name)
    : SubsysReco(name) {}

//____________________________________________________________________________..
SecondaryVertexTagger::~SecondaryVertexTagger() {}

//____________________________________________________________________________..
int SecondaryVertexTagger::Init(PHCompositeNode *topNode)
{
  if (m_truth_mode)
  {
    if (Verbosity() >= 1)
    {
      std::cout << __FILE__ << ": running in truth mode" << std::endl;
    }

    PHG4TruthInfoContainer *check_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (!check_truthInfo)
    {
      if (Verbosity() >= VERBOSITY_MORE)
       {
        std::cout << __FILE__ << ": Missing node G4TruthInfo" << std::endl;
        //exit(1);
      }
    }
  }

  if (m_save_output && Verbosity() >= 1)
  {
    std::cout << "Output nTuple: " << m_outfile_name << std::endl;
  }

  getField();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SecondaryVertexTagger::process_event(PHCompositeNode *topNode) {

  if (!m_truth_mode)
  { 
   SvtxTrackMap *check_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

    if (check_trackmap->size() == 0)
    {
      if (Verbosity() >= VERBOSITY_SOME)
      {
        std::cout << __FILE__ << ": Event skipped as there are no tracks and you're in reco mode" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  std::vector<KFParticle> tracks, vertices;

  tagVertices(topNode, tracks, vertices);

  if (tracks.size() != 0)
  {
    if (m_save_output && candidateCounter == 0)
    {
      m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
      initializeBranches();
    }

      candidateCounter += 1;

      if (m_save_output)
      {
        fillBranch(topNode, tracks, vertices, m_truth_mode);
      }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SecondaryVertexTagger::End(PHCompositeNode* /*topNode*/)
{
  if (m_save_output && candidateCounter != 0)
  {
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SecondaryVertexTagger::Reset(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
//void SecondaryVertexTagger::Print(const std::string &what) const {}

void SecondaryVertexTagger::getField() 
{
  // This sweeps the sPHENIX magnetic field map from some point radially then
  // grabs the first event that passes the selection
  m_magField = std::filesystem::exists(m_magField)
                   ? m_magField
                   : CDBInterface::instance()->getUrl(m_magField);

  if (Verbosity() > 0)
  {
    std::cout << __FILE__ << ": using fieldmap : " << m_magField << std::endl;
  }

  TFile *fin = new TFile(m_magField.c_str());
  TTree *fieldmap = (TTree *)fin->Get("fieldmap");

  float Bz = 0.;
  unsigned int r = 0.;
  float z = 0.;

  double arc = M_PI / 2;
  unsigned int n = 0;

  while (Bz == 0) 
  {
    if (n == 4)
    {
      ++r;
    }

    if (r == 3) // Dont go too far out radially
    {
      ++z;
    }

    n = n & 0x3U; // Constrains n from 0 to 3
    r = r & 0x2U;

    double x = r * std::cos(n * arc);
    double y = r * std::sin(n * arc);

    std::string sweep = "x == " + std::to_string(x) +
                        " && y == " + std::to_string(y) +
                        " && z == " + std::to_string(z);

    fieldmap->Draw(">>elist", sweep.c_str(), "entrylist");
    TEntryList *elist = (TEntryList *)gDirectory->Get("elist");
    if (elist->GetEntry(0))
    {
      TLeaf *fieldValue = fieldmap->GetLeaf("bz");
      fieldValue->GetBranch()->GetEntry(elist->GetEntry(0));
      Bz = fieldValue->GetValue();
    }

    ++n;

    if (r == 0) // No point in rescanning (0,0)
    {
      ++r;
      n = 0;
    }
  }
  // The actual unit of KFParticle is in kilo Gauss (kG), which is equivalent to
  // 0.1 T, instead of Tesla (T). The positive value indicates the B field is in
  // the +z direction
  Bz *= 10; // Factor of 10 to convert the B field unit from kG to T
  KFParticle::SetField((double)Bz);

  fieldmap->Delete();
  fin->Close();
}

