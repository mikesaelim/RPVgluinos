// Hadronizer.cc

// 8/27/13: This code implements either
//   1) the Pythia showering and hadronization of a generic sample of MC events provided as a LHE file, or
//   2) the Pythia event generation, showering, and hadronization for a set of Pythia initialization commands provided as a Pythia CMND file.
// The results are then output to a HepMC file, so that they can be passed on to Delphes or another program in the analysis chain.

// It is called with 3 cor 4 ommand-line arguments:
//    * (string)    path to the output HepMC file
//    * (string)    "LHE" if the input is a LHE file or "CMND" if the input is a CMND file
//   If the input is a LHE file,
//    * (string)    path to the input LHE file containing the MC sample from MadEvent
//   If the input is a CMND file, 
//    * (string)    path to the input CMND file containing the Pythia initialization commands, and
//    * (string)    path to the input SPC file containing the SLHA spectrum


// Output is to a HepMC file containing the fully hadronized event.

// This program must be compiled with the associated Makefile.





//////////////// Other program notes /////////////////

// This code is converted from /08-Snowmass2013/01-RPCNoDelphes/01-BoostedTopAnalysis.cc.

// Masses and energies in GeV.  Cross sections in pb.




#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <string>

#include "Pythia.h"
#include "HepMCInterface.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"




///////////////////////// Main function ////////////////////////////
int main(int argc, char *argv[]) {

  // Require 3 or 4 command-line arguments
  if(argc != 4 && argc != 5) {
    printf("!!!  ERROR: Wrong number of command-line arguments  !!!\n");
    printf("!!!     1st argument: path to the output HepMC file\n");
    printf("!!!     2nd argument: \"LHE\" if the input is a LHE file or \"CMND\" if the input is a CMND file\n");
    printf("!!!     3rd argument: path to the input LHE or CMND file\n");
    printf("!!!     4th argument: if it's a CMND file, path to the input SPC file\n");
    exit(1);
  }


  //************************** Initial setup, and retrieving command-line inputs **************************//
  // Run variables
  //    Inputs
  std::string const output_hepmc_filename = argv[1];      // Path to the output HepMC file
  std::string const input_type = argv[2];                 // "LHE" if the input is a LHE file or "CMND" if the input is a CMND file
  std::string const input_lhe_cmnd_filename = argv[3];    // Path to the input LHE or CMND file

  std::string input_spc_filename;                         // Path to the input SPC file
  if(input_type == "LHE") {
    input_spc_filename = "";
  }
  else if(input_type == "CMND") {
    input_spc_filename = argv[4];
  }
  else {
    printf("!!!  ERROR: Unknown input type, must be \"LHE\" or \"CMND\"  !!!\n");
    exit(1);
  }
  //******************************************************************************************************//


  printf("\n");
  printf("!!!!  Beginning hadronization on %s file:  %s  !!!!\n", input_type.c_str(), input_lhe_cmnd_filename.c_str());
  if(input_type == "CMND")  printf("!!!!     Spectrum file: %s  !!!!\n", input_spc_filename.c_str());
  printf("!!!!    Output to HepMC file:  %s  !!!!\n", output_hepmc_filename.c_str());
  printf("\n");

  
  //*********************************** Pythia initialization ***************************************//
  Pythia8::Pythia pythia;

  // Run parameters
  if(input_type == "LHE") {
    // The max number of events must be set by the user since the default is only 1K
    pythia.readString("Main:numberOfEvents = 1000000");  // Upper limit of 1M - program normally to be stopped by end of file
  }
  else {
    pythia.readString(std::string("SLHA:file = ") + input_spc_filename);
    pythia.readFile(input_lhe_cmnd_filename);
  }    
  pythia.readString("Main:timesAllowErrors = 10");

  // Use a new random number seed each time
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  // Extra settings
  pythia.readString("Init:showChangedSettings = off");      // Don't Show what settings were changed
  pythia.readString("Init:showChangedParticleData = off");  // Don't Show what particle data was changed
  pythia.readString("Next:numberCount = 1000");            // Show event count every 1000 events
  pythia.readString("Next:numberShowInfo = 0");            // Don't Show run info at beginning
  pythia.readString("Next:numberShowProcess = 0");         // Don't Show process record at beginning
  pythia.readString("Next:numberShowEvent = 0");           // Don't Show event record at beginning
  
  // For consistency, match PDF, renormalization and factorization scales with MadGraph defaults
  // Set PDF to CTEQ6L1
  pythia.readString("PDF:pSet = 8");
  // Set renormalization and factorization scales to the geometric mean of the transverse masses of the outgoing particles from the hard process
  pythia.readString("SigmaProcess:renormScale1 = 1");
  pythia.readString("SigmaProcess:renormScale2 = 2");
  pythia.readString("SigmaProcess:renormScale3 = 3");
  pythia.readString("SigmaProcess:factorScale1 = 1");
  pythia.readString("SigmaProcess:factorScale2 = 2");
  pythia.readString("SigmaProcess:factorScale3 = 3");
  
  // Finish initialization
  if(input_type == "LHE")  pythia.init(input_lhe_cmnd_filename);
  else                     pythia.init();
  //*************************************************************************************************//


  //*********************************** HepMC initialization ***************************************//
  // Initialize HepMC converter interface
  HepMC::I_Pythia8 hepmc_converter;
  // Initialize HepMC event output
  HepMC::IO_GenEvent hepmc_output(output_hepmc_filename, std::ios::out);
  //*************************************************************************************************//



  //////////////////////////////////////////////////////////////////////////////
  //                            Begin event loop                              //
  //////////////////////////////////////////////////////////////////////////////

  for(int i_event = 0; i_event < pythia.mode("Main:numberOfEvents"); ++i_event) {

    if(!pythia.next()) {  // Iterates Pythia to the next event
      if(pythia.info.atEndOfFile())  break;    // If at end of LHE file, stop

      printf("Warning: event %d failed\n", i_event);
      continue;
    }

    // Fill new HepMC GenEvent and pass it to the HepMC event output
    HepMC::GenEvent* hepmc_event = new HepMC::GenEvent();
    hepmc_converter.fill_next_event(pythia, hepmc_event);
    hepmc_output << hepmc_event;
    delete hepmc_event;
  }
  ///////////////////////////////////////////////////////////////////////////////
  //                           End event loop                                  //
  ///////////////////////////////////////////////////////////////////////////////

    
  // Done.                          
  return 0;
}






