#include "SBSSimEvent.h"
#include "TTree.h"
#include <iostream>

// -----------------------------------------------
// class SBSSimEvent: encapsulation of g4sbs_tree
//
/*
//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent() : g4sbs_tree()
{
  std::cout << "Initializing SBSSimEvent" << std::endl;
  RunID = EvtID = 0;
  //Weight = 1;
  Clear();
}
*/
//_____________________________________________________________________________
// SBSSimEvent::SBSSimEvent(TTree* tree, TString experiment) {
SBSSimEvent::SBSSimEvent(TTree *tree, Exp_t experiment) {
  std::cout << "Initializing SBSSimEvent" << std::endl;
  RunID = EvtID = 0;

  fExperiment = experiment;

  // Now we need to initialize the appropriate Tree structure based on experiment:
  // We should probably use an enum or something simple to make this less clunky than doing a string comparison each
  // time we open the file or load the event:
  // Tlad = new lad_tree_digitized(tree);
  Tgmn = new lad_tree_digitized(tree);
  // Tgep = new lad_tree_digitized(tree);
  // Tgenrp = new lad_tree_digitized(tree);

  // Weight = 1;
  Clear();
}
/*
//_____________________________________________________________________________
SBSSimEvent::SBSSimEvent(TTree* tree, std::vector<TString> det_list) : g4sbs_tree(tree, det_list)
{
  std::cout << "Initializing SBSSimEvent" << std::endl;
  cout << det_list.size() << endl;
  RunID = EvtID = 0;
  //Weight = 1;
  Clear();
}
*/
//_____________________________________________________________________________
void SBSSimEvent::Clear(const Option_t *opt) {
  // do nothing...
}

//_____________________________________________________________________________
void SBSSimEvent::Print(const Option_t *opt) const {
  // std::cout << RunID << " " << EvtID << " " << ev_sigma*ev_solang << std::endl;
}

//_____________________________________________________________________________
Int_t SBSSimEvent::GetEntry(Long64_t entry) {
  EvtID = entry;

  return Tgmn->GetEntry(entry);
  // return Tlad->GetEntry(entry);
}

//-----------------------------------------------------------------------------
ClassImp(SBSSimEvent)
