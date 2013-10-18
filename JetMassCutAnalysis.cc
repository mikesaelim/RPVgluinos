// JetMassCutAnalysis.cc

// 9/5/13: This code implements our CMS SSDL+b jets+MET search + high mass jet (HMJ) cut (arXiv:1302.2146) for the Snowmass 2013 proceedings, for 14 TeV and 300 fb^-1 and 3000 fb^-1.  This code handles the analysis steps after Delphes detector simulation.

// HepMC sample generation by Pythia is already handled by the Hadronizer program.
// The object selection cuts are handled by Delphes and specified in the accompanying Delphes card.
// Here event-level cuts, and the jet mass measurement and cut on the number of high-mass jets (> m(top)), are implemented.

// Notes from the earlier event generation part of the analysis chain:
//   For the signal process (gluino pair-production, with gluinos decaying to stop and top, and stops RPV-decaying to b and s), if m(stop) << m(gluino) - m(top), the program will generate on-shell stops (with stops 2-body decaying to bs), and if m(stop) >> m(gluino), the program will generate off-shell stops (with gluinos 3-body decaying to tbs).  As such, it is NOT advised to use this program for points where m(stop) ~ m(gluino), since that area of parameter space receives contributions from both on-shell stop and off-shell stop amplitudes.
//   Furthermore, I forced W's to decay leptonically for the signal process. This is the only way for the signal to generate two hard same-sign leptons.  So the SSDL cut partial efficiency will be multiplied by (0.3257)^2 to account for this.


// It is called with 5 command-line arguments:
//    * (string)  path to the output file
//    * (string)  process tag used to label this particular run/sample in the output (ex. "RPVgluino(800,300,50)", "ttbarW", etc.)
//    * (float)   LO cross section for this process, in pb
//    * (float)   K-factor for this process
//    * (string)  path to the input ROOT file containing the Delphes output


// It appends the output file with six sections:
//   The process tag and ROOT filename
//   The total LO production cross section, K-factor, total NLO production cross section, and total number of events generated for each point
//   For each signal region, the raw number of passing events and partial efficiency of the SSDL cut, of the # jets cut, of the # b jets cut, of the H_T cut, and of the MET cut
//   For each signal region, the total raw number of passing events and efficiency of all the CMS cuts, NLO cross section * efficiency, and resulting number of signal events for 300 and 3000 fb^-1.
//   For each signal region, the raw number of passing events and partial efficiency of the 1HMJ cut, and the total efficiency, NLO cross section * efficiency, and resulting number of signal events for 300 and 3000 fb^-1 for the CMS+1HMJ cuts.
//   For each signal region, the raw number of passing events and partial efficiency of the 2HMJ cut, and the total efficiency, NLO cross section * efficiency, and resulting number of signal events for 300 and 3000 fb^-1 for the CMS+2HMJ cuts.


// This program must be compiled with the associated Makefile.





//////////////// Other program notes /////////////////

// This code is converted from /08-Snowmass2013/02-RPVNoDelphes/01-JetMassCutAnalysis.cc and /07-StealthStopSSDL/02-CMS-SUS-13-013/01-CMSrecast.cc, with Delphes ROOT TTree access learned from /Delphes-3.0.10/examples/Example1.cpp and Example1.C.
// Lepton identification efficiencies were calculated in /04-.../05-.../02-ReplicateLM9.cc: 73% for electrons and 84% for muons.

// Masses and energies in GeV.  Cross sections in pb.





#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <string>
#include <utility>

#include "fastjet/PseudoJet.hh"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"




//////////////////////////// Helper struct and function declarations /////////////////////////////
// This unnamed namespace conveniently encapsulates helper structs and free helper functions.
// Functions are implemented below the main function.
namespace {
  // Struct for storing signal region definitions and results
  //   Reasons for using a struct instead of a class:
  //     1) It might be nice to protect the SR settings from being changed after initialization.
  //        However, using const to do that means that the default copy/assignment operators no longer work, and those operators are needed
  //          by std::vector::push_back.
  //        But constructing a proper class for SignalRegion is complete overkill if it's just for this program.
  //     2) Initializing with ={...} is the quickest way to load all the signal regions at the beginning of the analysis
  struct SignalRegion {
    unsigned int num_jets;     // Minimum number of jets passing kinematic cuts
    unsigned int num_bjets;    // Minimum number of tagged b jets
    int ssdl_charge;         // Given SSDL, 1 to require ++, -1 to require --, and 0 to allow for both ++/--
    double ht_cut;             // Minimum H_T
    double met_cut;            // Minimum MET

    unsigned int num_events_ssdl_cut;            // Raw number of events passing the SSDL cuts
    double part_eff_ssdl_cut;                    // Partial efficiency of the SSDL cuts
    unsigned int num_events_jet_cut;             // Raw number of events passing the SSDL + jets cuts
    double part_eff_jet_cut;                     // Partial efficiency of the jets cut
    unsigned int num_events_bjet_cut;            // Raw number of events passing the SSDL + jets + b jets cuts
    double part_eff_bjet_cut;                    // Partial efficiency of the b jets cut
    unsigned int num_events_ht_cut;              // Raw number of events passing the SSDL + jets + b jets + H_T cuts
    double part_eff_ht_cut;                      // Partial efficiency of the H_T cut
    unsigned int num_events_met_cut;             // Raw number of events passing the SSDL + jets + b jets + H_T + MET cuts
    double part_eff_met_cut;                     // Partial efficiency of the MET cut

    unsigned int num_events_cms;                 // Raw number of events passing all the CMS cuts
    double total_eff_cms;                        // Total efficiency of the CMS cuts
    double nlo_cross_section_passed_cms;         // NLO cross section * total efficiency of the CMS cuts
    double num_events_for_int_lumi1_cms;         // Number of events expected for integrated luminosity int_lumi1, for the CMS cuts
    double num_events_for_int_lumi2_cms;         // Number of events expected for integrated luminosity int_lumi2, for the CMS cuts

    unsigned int num_events_1hmj;          // Raw number of events passing the CMS+1HMJ cuts
    double part_eff_1hmj;                  // Partial efficiency of the 1HMJ cut
    double total_eff_1hmj;                 // Total efficiency of the CMS+1HMJ cuts
    double nlo_cross_section_passed_1hmj;  // NLO cross section * total efficiency of the CMS+1HMJ cuts 
    double num_events_for_int_lumi1_1hmj;  // Number of events expected for integrated luminosity int_lumi1, for the CMS+1HMJ cuts
    double num_events_for_int_lumi2_1hmj;  // Number of events expected for integrated luminosity int_lumi2, for the CMS+1HMJ cuts

    unsigned int num_events_2hmj;          // Raw number of events passing the CMS+2HMJ cuts
    double part_eff_2hmj;                  // Partial efficiency of the 2HMJ cut
    double total_eff_2hmj;                 // Total efficiency of the CMS+2HMJ cuts
    double nlo_cross_section_passed_2hmj;  // NLO cross section * total efficiency of the CMS+2HMJ cuts
    double num_events_for_int_lumi1_2hmj;  // Number of events expected for integrated luminosity int_lumi1, for the CMS+2HMJ cuts
    double num_events_for_int_lumi2_2hmj;  // Number of events expected for integrated luminosity int_lumi2, for the CMS+2HMJ cuts
  };

  //////// Helper functions ////////
  std::vector<SignalRegion> CreateAssay();                     // Generates all the signal regions in the analysis, returns completed assay

  //   Small utility functions
  bool IsSamePseudojet(fastjet::PseudoJet pseudojet1, fastjet::PseudoJet pseudojet2);      // returns true if the fastjet::PseudoJets are equal
  double PairInvariantMass(fastjet::PseudoJet pseudojet1, fastjet::PseudoJet pseudojet2);  // returns the pair invariant mass of two fastjet::PseudoJets
  int GetSign(int num);                                                                    // returns the sign (-1, 0, +1) of the integer
  double RandNum();                                                                        // returns a random number between 0 and 1
  
  //   SSDL-oriented functions
  bool ApplyDileptonTrigger(int id1, int id2);                                             // Apply dilepton trigger efficiencies
  bool ApplyThirdLeptonVeto(std::pair<fastjet::PseudoJet, fastjet::PseudoJet> ssdl_pair, 
			    std::vector<fastjet::PseudoJet> loose_leptons,
			    double third_lepton_pt_cut,
			    double pair_mass_min, double pair_mass_max);                   // Apply vetoes on 3rd leptons forming resonances with SSDL leptons
                                                                                           // Note: returns FALSE if the event is vetoed!

  bool TagHMJet(double JetMass);                                                           // Tag high mass jets (mass above m(top))   
}
//////////////////////////////////////////////////////////////////////////////////////////////////






int main(int argc, char *argv[]) {

  // Require 5 command-line arguments
  if(argc != 6) {
    printf("Wrong number of command-line arguments!\n");
    exit(1);
  }
  

  //************************** Initial setup, and retrieving command-line inputs **************************//
  // Run variables
  //    Inputs
  std::string const output_filename = argv[1];             // Path to the output file
  std::string const process_tag = argv[2];                 // Label for this particular run/sample in the output
  double const lo_cross_section = atof(argv[3]);           // LO cross section in pb
  double const k_factor = atof(argv[4]);                   // K factor
  std::string const root_filename = argv[5];               // Path to the ROOT file
  //     Outputs
  double nlo_cross_section = lo_cross_section * k_factor;  // NLO cross section in pb
  unsigned int num_events = 0;                             // Raw number of events processed at this point
  std::vector<SignalRegion> assay = CreateAssay();         // Vector of signal regions - terminology borrowed from biology

  // The RPV signal will need special handling for its SSDL cut efficiency
  bool is_signal;
  if(process_tag.find("RPVgluinos") != std::string::npos)  is_signal = true;
  else                                                     is_signal = false;

  // Fixed, hard-coded parameters
  double const int_lumi1 = 300.;              // Integrated luminosities in fb^-1
  double const int_lumi2 = 3000.;

  // Seed random number generator with the current time
  srand( time(NULL) );
  //********************************************************************************************************//


  printf("\n");
  printf("!!!!  Beginning analysis on %s   !!!!\n", process_tag.c_str());
  printf("!!!!      Input from %s       !!!!\n", root_filename.c_str());
  printf("!!!!      Output to %s       !!!!\n", output_filename.c_str());
  printf("\n");  



  //**************************************** Prepare input ROOT TTree *******************************************//
  // Create ExRootTreeReader from ROOT file
  TChain input_chain("Delphes");
  input_chain.Add(root_filename.c_str());
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&input_chain);
  
  // Get pointers to branches
  //    Main object collections
  TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
  TClonesArray *branch_electron = tree_reader->UseBranch("Electron");
  TClonesArray *branch_muon = tree_reader->UseBranch("Muon");
  TClonesArray *branch_ht = tree_reader->UseBranch("ScalarHT");
  TClonesArray *branch_met = tree_reader->UseBranch("MissingET");
  //    Loose object collections for event vetoes
  TClonesArray *branch_loose_electron = tree_reader->UseBranch("looseElectron");
  TClonesArray *branch_loose_muon = tree_reader->UseBranch("looseMuon"); 
  //*************************************************************************************************************//
  


  //////////////////////////////////////////////////////////////////////////////
  //                            Begin event loop                              //
  //////////////////////////////////////////////////////////////////////////////
  
  for(int i_event = 0; i_event < tree_reader->GetEntries(); ++i_event) {
    
    // Increment the number of events generated
    ++num_events;

    // Load branches with data from the event
    tree_reader->ReadEntry(i_event);
    

    //************************************* Extract the objects from the branches **************************************//
    // Here I also double-check the pT and eta cuts for each object
    // In order to process the leptons as I had before, I convert the electrons and muons into fastjet::PseudoJet objects
    // The user_index of the fastjet::PseudoJet objects is the PDG ID of the lepton
    std::vector<fastjet::PseudoJet> leptons;
    for(int i_electron = 0; i_electron < branch_electron->GetEntries(); ++i_electron) {
      Electron *delphes_electron = dynamic_cast<Electron*>(branch_electron->At(i_electron));
      if(!(delphes_electron->PT >= 20.))  continue;
      if(!(abs(delphes_electron->Eta) <= 1.442 || (abs(delphes_electron->Eta) >= 1.566 && abs(delphes_electron->Eta) <= 2.4)))  continue;

      TLorentzVector vec = delphes_electron->P4();
      fastjet::PseudoJet pj_electron(vec.Px(), vec.Py(), vec.Pz(), vec.E());
      pj_electron.set_user_index(delphes_electron->Charge * 11);

      leptons.push_back(pj_electron);
    }
    for(int i_muon = 0; i_muon < branch_muon->GetEntries(); ++i_muon) {
      Muon *delphes_muon = dynamic_cast<Muon*>(branch_muon->At(i_muon));
      if(!(delphes_muon->PT >= 20.))  continue;
      if(!(abs(delphes_muon->Eta) <= 2.4))  continue;

      TLorentzVector vec = delphes_muon->P4();
      fastjet::PseudoJet pj_muon(vec.Px(), vec.Py(), vec.Pz(), vec.E());
      pj_muon.set_user_index(delphes_muon->Charge * 13);

      leptons.push_back(pj_muon);
    }

    std::vector<Jet*> jets;
    std::vector<Jet*> bjets;
    std::vector<Jet*> hmjets;
    for(int i_jet = 0; i_jet < branch_jet->GetEntries(); ++i_jet) {
      Jet *delphes_jet = dynamic_cast<Jet*>(branch_jet->At(i_jet));
      if(!(delphes_jet->PT >= 40.))  continue;
      if(!(abs(delphes_jet->Eta) <= 2.4))  continue;

      jets.push_back(delphes_jet);
      if(delphes_jet->BTag == 1)  bjets.push_back(delphes_jet);
      if(TagHMJet(delphes_jet->Mass))  hmjets.push_back(delphes_jet);
    }

    double ht = dynamic_cast<ScalarHT*>(branch_ht->At(0))->HT;
    double met = dynamic_cast<MissingET*>(branch_met->At(0))->MET;

    std::vector<fastjet::PseudoJet> loose_leptons;
    for(int i_electron = 0; i_electron < branch_loose_electron->GetEntries(); ++i_electron) {
      Electron *delphes_electron = dynamic_cast<Electron*>(branch_loose_electron->At(i_electron));
      if(!(delphes_electron->PT >= 20.))  continue;
      if(!(abs(delphes_electron->Eta) <= 1.442 || (abs(delphes_electron->Eta) >= 1.566 && abs(delphes_electron->Eta) <= 2.4)))  continue;

      TLorentzVector vec = delphes_electron->P4();
      fastjet::PseudoJet pj_electron(vec.Px(), vec.Py(), vec.Pz(), vec.E());
      pj_electron.set_user_index(delphes_electron->Charge * 11);

      loose_leptons.push_back(pj_electron);
    }
    for(int i_muon = 0; i_muon < branch_loose_muon->GetEntries(); ++i_muon) {
      Muon *delphes_muon = dynamic_cast<Muon*>(branch_loose_muon->At(i_muon));
      if(!(delphes_muon->PT >= 20.))  continue;
      if(!(abs(delphes_muon->Eta) <= 2.4))  continue;

      TLorentzVector vec = delphes_muon->P4();
      fastjet::PseudoJet pj_muon(vec.Px(), vec.Py(), vec.Pz(), vec.E());
      pj_muon.set_user_index(delphes_muon->Charge * 13);
      
      loose_leptons.push_back(pj_muon);
    }
    //************************************************************************************************************//
    
    

    //************************************** SSDL cut analysis ****************************************//
    // Apply SSDL cut
    // Sort selected leptons by pT, descending 
    //   (simple bubble sort - only a few leptons anyway)
    //   (can use swap function - fastjet::PseudoJet is ok with default copy and assignment constructors)
    bool did_swap = true;
    while(did_swap) {
      did_swap = false;
      for(unsigned int i_lepton = 0; i_lepton + 1 < leptons.size(); ++i_lepton) {
	if(leptons[i_lepton].pt() < leptons[i_lepton+1].pt()) {
	  std::swap(leptons[i_lepton], leptons[i_lepton+1]);
	  did_swap = true;
	}
      }
    }
    // Look for the pair of SS leptons with the highest leading pT, and call that our dilepton candidates.  We also require the invariant mass > 8 GeV.
    bool ssdl_found = false;
    std::pair<fastjet::PseudoJet, fastjet::PseudoJet> ssdl_pair;
    for(std::vector<fastjet::PseudoJet>::const_iterator i_lepton1 = leptons.begin(); i_lepton1 < leptons.end(); ++i_lepton1) {
      for(std::vector<fastjet::PseudoJet>::const_iterator i_lepton2 = i_lepton1 + 1; i_lepton2 < leptons.end(); ++i_lepton2) {
	if(GetSign(i_lepton1->user_index()) == GetSign(i_lepton2->user_index()) &&
	   PairInvariantMass(*i_lepton1, *i_lepton2) > 8.) {
	  ssdl_found = true;
	  ssdl_pair.first = *i_lepton1;
	  ssdl_pair.second = *i_lepton2;
	  break;
	}
      }
      if(ssdl_found)  break;
    }
    if(!ssdl_found)  continue;
    
    // Apply dilepton trigger efficiency
    if(!ApplyDileptonTrigger(ssdl_pair.first.user_index(), ssdl_pair.second.user_index()))  continue;
    
    // Veto events where a third loose lepton (w/ pT > 10 GeV) forms a possible Z pair with one of the dilepton candidates
    if(!ApplyThirdLeptonVeto(ssdl_pair, loose_leptons, 10., 76., 106.))  continue;
    
    // Veto events where a third loose lepton (w/ pT > 5 GeV) forms a possible gamma* or low mass neutral bound state pair with one of the dilepton candidates
    if(!ApplyThirdLeptonVeto(ssdl_pair, loose_leptons, 5., 0., 12.))  continue;
    //**********************************************************************************************//
    
      
      
    //**************************** Apply event-level cuts for each signal region ******************************//
    // SSDL cut (except for pair charge sign) already applied in the leptonic analysis
    for(std::vector<SignalRegion>::iterator i_sr = assay.begin(); i_sr < assay.end(); ++i_sr) {
      if(i_sr->ssdl_charge == 0)  ++i_sr->num_events_ssdl_cut;
      else {
	if(GetSign(ssdl_pair.first.user_index()) == i_sr->ssdl_charge)  ++i_sr->num_events_ssdl_cut;
	else  continue;
      }

      if(jets.size() >= i_sr->num_jets)  ++i_sr->num_events_jet_cut;
      else  continue;

      if(bjets.size() >= i_sr->num_bjets)  ++i_sr->num_events_bjet_cut;
      else  continue;

      if(ht >= i_sr->ht_cut)  ++i_sr->num_events_ht_cut;
      else  continue;

      if(met >= i_sr->met_cut)  ++i_sr->num_events_met_cut;
      else  continue;

      ++i_sr->num_events_cms;
      
      if(hmjets.size() >= 1)  ++i_sr->num_events_1hmj;
      else  continue;

      if(hmjets.size() >= 2)  ++i_sr->num_events_2hmj;
      else  continue;
    }
    //********************************************************************************************//
  }
  ///////////////////////////////////////////////////////////////////////////////
  //                           End event loop                                  //
  ///////////////////////////////////////////////////////////////////////////////

  
  
  //*********************************** Begin run analysis ************************************//
  // Parse the results for each signal region - casting both numerator and denominator so that 0 in denominator is properly handled
  for(std::vector<SignalRegion>::iterator i_sr = assay.begin(); i_sr < assay.end(); ++i_sr) {
    i_sr->part_eff_ssdl_cut = static_cast<double>(i_sr->num_events_ssdl_cut) / static_cast<double>(num_events);
    if(is_signal)  i_sr->part_eff_ssdl_cut *= pow(0.3257, 2);   // Correction for forcing W's to decay leptonically
    i_sr->part_eff_jet_cut = static_cast<double>(i_sr->num_events_jet_cut) / static_cast<double>(i_sr->num_events_ssdl_cut);
    i_sr->part_eff_bjet_cut = static_cast<double>(i_sr->num_events_bjet_cut) / static_cast<double>(i_sr->num_events_jet_cut);
    i_sr->part_eff_ht_cut = static_cast<double>(i_sr->num_events_ht_cut) / static_cast<double>(i_sr->num_events_bjet_cut);
    i_sr->part_eff_met_cut = static_cast<double>(i_sr->num_events_met_cut) / static_cast<double>(i_sr->num_events_ht_cut);

    i_sr->total_eff_cms = static_cast<double>(i_sr->num_events_cms) / static_cast<double>(num_events);
    if(is_signal)  i_sr->total_eff_cms *= pow(0.3257, 2);   // Correction for forcing W's to decay leptonically
    i_sr->nlo_cross_section_passed_cms = nlo_cross_section * i_sr->total_eff_cms;
    i_sr->num_events_for_int_lumi1_cms = i_sr->nlo_cross_section_passed_cms * int_lumi1 * 1.E+03;
    i_sr->num_events_for_int_lumi2_cms = i_sr->nlo_cross_section_passed_cms * int_lumi2 * 1.E+03;

    i_sr->part_eff_1hmj = static_cast<double>(i_sr->num_events_1hmj) / static_cast<double>(i_sr->num_events_cms);
    i_sr->total_eff_1hmj = static_cast<double>(i_sr->num_events_1hmj) / static_cast<double>(num_events);
    if(is_signal)  i_sr->total_eff_1hmj *= pow(0.3257, 2);   // Correction for forcing W's to decay leptonically
    i_sr->nlo_cross_section_passed_1hmj = nlo_cross_section * i_sr->total_eff_1hmj;
    i_sr->num_events_for_int_lumi1_1hmj = i_sr->nlo_cross_section_passed_1hmj * int_lumi1 * 1.E+03;
    i_sr->num_events_for_int_lumi2_1hmj = i_sr->nlo_cross_section_passed_1hmj * int_lumi2 * 1.E+03;

    i_sr->part_eff_2hmj = static_cast<double>(i_sr->num_events_2hmj) / static_cast<double>(i_sr->num_events_1hmj);
    i_sr->total_eff_2hmj = static_cast<double>(i_sr->num_events_2hmj) / static_cast<double>(num_events);
    if(is_signal)  i_sr->total_eff_2hmj *= pow(0.3257, 2);   // Correction for forcing W's to decay leptonically
    i_sr->nlo_cross_section_passed_2hmj = nlo_cross_section * i_sr->total_eff_2hmj;
    i_sr->num_events_for_int_lumi1_2hmj = i_sr->nlo_cross_section_passed_2hmj * int_lumi1 * 1.E+03;
    i_sr->num_events_for_int_lumi2_2hmj = i_sr->nlo_cross_section_passed_2hmj * int_lumi2 * 1.E+03;
  }
  //*******************************************************************************************//
  


  //************************************* Print the results ***************************************//
  printf("\n");
  printf("  Printing results to file...\n");
  
  FILE* output_file = fopen(output_filename.c_str(), "a");

  fprintf(output_file, "  process                                    ROOT file\n");
  fprintf(output_file, " %s             %s\n", process_tag.c_str(), root_filename.c_str());
  fprintf(output_file, "  sigmaLO (pb)        K-factor         sigmaNLO (pb)       N\n");
  fprintf(output_file, "%- 9.8E    %- 9.8E    %- 9.8E    %u\n", lo_cross_section, k_factor, nlo_cross_section, num_events);

  fprintf(output_file, " SR#     N(SSDL cut)     eff(SSDL cut)     N(jet cut)     eff(jet cut)     N(bjet cut)     eff(bjet cut)     N(H_T cut)     eff(H_T cut)     N(MET cut)     eff(MET cut)\n");
  for(unsigned int i_sr = 0; i_sr < assay.size(); ++i_sr) {
    fprintf(output_file, "%3u       %6u       %- 9.8E      %6u      %- 9.8E     %6u        %- 9.8E    %6u       %- 9.8E     %6u      %- 9.8E\n", i_sr, assay[i_sr].num_events_ssdl_cut, assay[i_sr].part_eff_ssdl_cut, assay[i_sr].num_events_jet_cut, assay[i_sr].part_eff_jet_cut, assay[i_sr].num_events_bjet_cut, assay[i_sr].part_eff_bjet_cut, assay[i_sr].num_events_ht_cut, assay[i_sr].part_eff_ht_cut, assay[i_sr].num_events_met_cut, assay[i_sr].part_eff_met_cut);
  }
  
  fprintf(output_file, " SR#     N(CMS cuts)     eff(CMS cuts)    sigmaNLO*eff(CMS cuts) (pb)    # events in %g fb^-1     # events in %g fb^-1\n", int_lumi1, int_lumi2);
  for(unsigned int i_sr = 0; i_sr < assay.size(); ++i_sr) {
    fprintf(output_file, "%3u       %6u       %- 9.8E        %- 9.8E              %- 9.8E            %- 9.8E\n", i_sr, assay[i_sr].num_events_cms, assay[i_sr].total_eff_cms, assay[i_sr].nlo_cross_section_passed_cms, assay[i_sr].num_events_for_int_lumi1_cms, assay[i_sr].num_events_for_int_lumi2_cms);
  }

  fprintf(output_file, " SR#     N(1HMJ cut)     eff(1HMJ cut)     eff(CMS+1HMJ)    sigmaNLO*eff(CMS+1HMJ) (pb)    # events in %g fb^-1     # events in %g fb^-1\n", int_lumi1, int_lumi2);
  for(unsigned int i_sr = 0; i_sr < assay.size(); ++i_sr) {
    fprintf(output_file, "%3u       %6u       %- 9.8E    %- 9.8E       %- 9.8E              %- 9.8E            %- 9.8E\n", i_sr, assay[i_sr].num_events_1hmj, assay[i_sr].part_eff_1hmj, assay[i_sr].total_eff_1hmj, assay[i_sr].nlo_cross_section_passed_1hmj, assay[i_sr].num_events_for_int_lumi1_1hmj, assay[i_sr].num_events_for_int_lumi2_1hmj);
  }

  fprintf(output_file, " SR#     N(2HMJ cut)     eff(2HMJ cut)     eff(CMS+2HMJ)    sigmaNLO*eff(CMS+2HMJ) (pb)    # events in %g fb^-1     # events in %g fb^-1\n", int_lumi1, int_lumi2);
  for(unsigned int i_sr = 0; i_sr < assay.size(); ++i_sr) {
    fprintf(output_file, "%3u       %6u       %- 9.8E    %- 9.8E       %- 9.8E              %- 9.8E            %- 9.8E\n", i_sr, assay[i_sr].num_events_2hmj, assay[i_sr].part_eff_2hmj, assay[i_sr].total_eff_2hmj, assay[i_sr].nlo_cross_section_passed_2hmj, assay[i_sr].num_events_for_int_lumi1_2hmj, assay[i_sr].num_events_for_int_lumi2_2hmj);
  }
  
  fprintf(output_file, "\n");
  
  fclose(output_file);
  //******************************************************************************************************//

  
  // Done.                          
  return 0;
}










///////////////////////// Helper functions, enclosed in unnamed namespace ///////////////////////////////
namespace {

  std::vector<SignalRegion> CreateAssay() {
    // Generates all the signal regions in the analysis, returns completed assay
    // Signal regions: # jets, # b jets, SSDL charge, H_T, MET, followed by results

    std::vector<SignalRegion> assay;

    SignalRegion sr0 = {2, 2, 0, 80.,  0.,   0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr1 = {2, 2, 0, 80.,  30.,  0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr2 = {2, 2, 1, 80.,  30.,  0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr3 = {4, 2, 0, 200., 120., 0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr4 = {4, 2, 0, 200., 50.,  0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr5 = {4, 2, 0, 320., 50.,  0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr6 = {4, 2, 0, 320., 120., 0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr7 = {3, 3, 0, 200., 50.,  0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    SignalRegion sr8 = {4, 2, 0, 320., 0.,   0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,0,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1.,0,-1.,-1.,-1.,-1.,-1};
    assay.push_back(sr0);
    assay.push_back(sr1);
    assay.push_back(sr2);
    assay.push_back(sr3);
    assay.push_back(sr4);
    assay.push_back(sr5);
    assay.push_back(sr6);
    assay.push_back(sr7);
    assay.push_back(sr8);
    
    return assay;    
  }



  bool IsSamePseudojet(fastjet::PseudoJet pseudojet1, fastjet::PseudoJet pseudojet2) {
    // returns true if the fastjet::PseudoJets are equal
    
    return ( pseudojet1.px() == pseudojet2.px() && 
	     pseudojet1.py() == pseudojet2.py() && 
	     pseudojet1.pz() == pseudojet2.pz() && 
	     pseudojet1.e() == pseudojet2.e() && 
	     pseudojet1.user_index() == pseudojet2.user_index() );
  }
  
  
  
  double PairInvariantMass(fastjet::PseudoJet pseudojet1, fastjet::PseudoJet pseudojet2) {
    // returns the pair invariant mass of two fastjet::PseudoJets
    
    fastjet::PseudoJet sum(pseudojet1);
    sum += pseudojet2;
    
    return sum.m();
  }
  
  
  
  int GetSign(int num) {
    // returns the sign (-1, 0, +1) of the integer
    if(num == 0)  return 0;
    
    return num / abs(num);
  }
  


  double RandNum() {
    // returns a random number between 0 and 1
    return static_cast<double>(rand()) / RAND_MAX;
  }
  
  
  
  bool ApplyDileptonTrigger(int id1, int id2) {
    // Apply trigger efficiencies for the dilepton trigger
    // Returns true if the event passes the dilepton trigger
    // Assumes that the event has passed the SSDL cut, so that the two leptons already pass the dilepton trigger's pT cuts
    
    if((abs(id1) != 11 && abs(id1) != 13) || (abs(id2) != 11 && abs(id2) != 13)) {
      printf("!!! ERROR: Particle IDs %d and %d passed to ApplyDileptonTrigger()!\n", id1, id2);
      exit(1);
    }
    
    double eff;
    if     (abs(id1) + abs(id2) == 22)  eff = 0.96;
    else if(abs(id1) + abs(id2) == 24)  eff = 0.93;
    else if(abs(id1) + abs(id2) == 26)  eff = 0.88;
    else                                eff = 0;
    
    return (RandNum() <= eff);
  }
  
  
  
  bool ApplyThirdLeptonVeto(std::pair<fastjet::PseudoJet, fastjet::PseudoJet> ssdl_pair, std::vector<fastjet::PseudoJet> loose_leptons, 
			    double third_lepton_pt_cut, double pair_mass_min, double pair_mass_max) {
    // Apply vetoes on 3rd leptons forming resonances with SSDL leptons
    // Returns FALSE if the event is vetoed
    
    for(std::vector<fastjet::PseudoJet>::const_iterator i_third_lepton = loose_leptons.begin(); i_third_lepton < loose_leptons.end(); ++i_third_lepton) {
      if(IsSamePseudojet(*i_third_lepton, ssdl_pair.first) || IsSamePseudojet(*i_third_lepton, ssdl_pair.second))  continue;
      
      // Apply third lepton pT cut
      if(i_third_lepton->pt() < third_lepton_pt_cut)  continue;
      
      // If the third lepton forms an opposite-sign, same-flavor pair with one of the dilepton candidates, with invariant mass between the limits, veto event.
      if(i_third_lepton->user_index() == -1 * ssdl_pair.first.user_index()) {
	if(PairInvariantMass(*i_third_lepton, ssdl_pair.first) > pair_mass_min && PairInvariantMass(*i_third_lepton, ssdl_pair.first) < pair_mass_max) {
	  return false;
	}
      }
      if(i_third_lepton->user_index() == -1 * ssdl_pair.second.user_index()) {
	if(PairInvariantMass(*i_third_lepton, ssdl_pair.second) > pair_mass_min && PairInvariantMass(*i_third_lepton, ssdl_pair.second) < pair_mass_max) {
	  return false;
	}
      }
    }
    
    return true;
  }



  bool TagHMJet(double JetMass) {
    // Tag high mass jets (mass above m(top))
    // Returns true if JetMass is above the top mass (PDG 2012 central value)
    
    return (JetMass >= 173.5);
  }
  
    

}





