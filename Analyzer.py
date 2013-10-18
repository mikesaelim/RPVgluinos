#!/usr/local/bin/python

""" Analyzer.py - Runs the full analysis chain for the RPV gluinos search for Snowmass 2013.

This script ties all the steps of the analysis together, taking in the scan parameters and executing the scan over the RPV model parameter space.
* Analysis: SSDL+b jets+MET search at 8 TeV and 10.5 fb^-1 (CMS SUS-12-017) with R = 1.0 and a cut on high mass jets (arXiv:1302.2146).
* Collider inputs: only 14 TeV at 300 and 3000 fb^-1, for Snowmass 2013.

* Signal event generation, and both signal and background hadronization, is handled with PYTHIA 8 by the model-specific GenerateInputs program and the general-purpose Hadronizer program.
* Detector simulation and object reconstruction is handled by Delphes 3 and specified by Delphes/delphes_card_CMS_SUS-12-017_R1mod.tcl
* Event-level cuts and analysis is handled by JetMassCutAnalysis.

For each point in the signal parameter space, the analysis chain is:
[param space point] -(GenerateInputs)-> [.spc,.cmnd] -(Hadronizer)-> [.hepmc] -(Delphes)-> [.root] -(JetMassCutAnalysis)-> [scan results .dat]
For each background sample, the analysis chain is:
[background info] -(Hadronizer)-> [.hepmc] -(Delphes)-> [.root] -(JetMassCutAnalysis)-> [scan results .dat]

It should be noted that the .root file that Delphes outputs does not contain the LO cross section, so the LO cross section must be extracted from the .hepmc file and delivered to JetMassCutAnalysis separately.  Afterwards, the .hepmc file is deleted because keeping it would use up too much hard drive space.


########## How this script, specifically, works ##########

Scan parameters are hard-coded near the beginning of the script.  The user has two choices for the run:
   * scan over specific background samples
   * grid scan over the signal parameter space with some user-defined limits

It takes in signal K-factors from a separate file.  If a K-factor isn't found, the point is skipped.
If a sample label is already present in the output file, either because the sample was already analyzed or because the label was already used by another run, the point is skipped.



--Written by Mike Saelim (mjs496@cornell.edu)

"""


import os
import sys
import time


################################################### BEGIN USER SETTINGS #########################################################

### scan_description is used to describe the run in the log
scan_description = "CM energy = 14 TeV @ 300 and 300 fb^-1; 10k events/point; m(gluino) = 1950; m(stop) = [1200,1600,100]"
#10k events/point; m(gluino) = [1500,2400,300]; m(stop) = [200,m(gluino)-300,400]"
#100k event backgrounds ttbarW_1 and ttbarZ_1"

### Filepaths - specify relative filepaths relative to the script's directory!
# Input generator executable
geninputs_exec = "./GenerateInputs"
# Hadronizer executable
hadronizer_exec = "./Hadronizer"
# Delphes executable (HepMC interface)
delphes_exec = "/Users/donerkebab/Desktop/work/software/Delphes-3.0.10/DelphesHepMC"
# Delphes card
delphes_card_file = "Delphes/delphes_card_CMS_SUS-12-017_R1mod.tcl"
# Final analysis executable
analysis_exec = "./JetMassCutAnalysis"
# Storing results
results_filename = "Results/JetMassCut_signals_14TeV.dat"   #"Results/JetMassCut_backgrounds_14TeV.dat"
hepmc_file_folder = "Results/HepMC"
root_file_folder = "Results/ROOT"
logs_folder = "logs"

### General fixed parameters
# Program only attempts to run each step of the analysis on a given sample this many times
max_attempts_per_sample = 3



### Scan parameters

# Type of scan to do
scan_mode = "signal"    # "signal" for a signal grid scan, "background" for a specific set of backgrounds


# If we're scanning through a specific set of backgrounds, specify them here as (sample_label, K-factor, path to gzipped LHE file)
# ttbarW: 1.236 (arXiv:1204.5678)
# ttbarZ: 1.387 (arXiv:1111.0610)
if scan_mode == "background":
    backgrounds_list = [
	("ttbarW_SSDL_14TeV_100k_1", 1.236, "Backgrounds/ttbarW_SSDL_14TeV_100k_1.lhe.gz"),
	("ttbarW_SSDL_14TeV_100k_2", 1.236, "Backgrounds/ttbarW_SSDL_14TeV_100k_2.lhe.gz"),
	("ttbarW_SSDL_14TeV_100k_3", 1.236, "Backgrounds/ttbarW_SSDL_14TeV_100k_3.lhe.gz"),
	("ttbarW_SSDL_14TeV_100k_4", 1.236, "Backgrounds/ttbarW_SSDL_14TeV_100k_4.lhe.gz"),
	("ttbarW_SSDL_14TeV_100k_5", 1.236, "Backgrounds/ttbarW_SSDL_14TeV_100k_5.lhe.gz"),
	("ttbarZ_SSDL_14TeV_100k_1", 1.387, "Backgrounds/ttbarZ_SSDL_14TeV_100k_1.lhe.gz"),
	("ttbarZ_SSDL_14TeV_100k_2", 1.387, "Backgrounds/ttbarZ_SSDL_14TeV_100k_2.lhe.gz"),
	("ttbarZ_SSDL_14TeV_100k_3", 1.387, "Backgrounds/ttbarZ_SSDL_14TeV_100k_3.lhe.gz"),
	("ttbarZ_SSDL_14TeV_100k_4", 1.387, "Backgrounds/ttbarZ_SSDL_14TeV_100k_4.lhe.gz"),
	("ttbarZ_SSDL_14TeV_100k_5", 1.387, "Backgrounds/ttbarZ_SSDL_14TeV_100k_5.lhe.gz")
        ]


# If we're doing a grid scan over the signal parameter space, specify...
if scan_mode == "signal":
    # Number of events to generate for each point
    num_events = 10000
    # K-factors
    kfactors_filename = "Kfactors_GPP_14TeV.dat"
    
    # Function for generating sample labels automatically
    def get_sample_label(gluino_mass, stop_mass):
        return "RPVgluinos_" + str(gluino_mass) + "_" + str(stop_mass)

    # Loops for generating the grid
    #   The terminating number in range() needs to be greater than the actual stop value
    signals_list = []
    for gluino_mass in range(1950, 1950 + 1, 50):
	for stop_mass in range(1200, 1600 + 1, 100):
	    signals_list.append((gluino_mass, stop_mass))


###################################################### END USER SETTINGS ###########################################################


### Change the working directory to the script's location
os.chdir(os.path.dirname(sys.argv[0]))

### Error count
error_count = 0

### Open log file and set up dual output
log_file_index = 1
log_filename = ""
while True:
    if os.path.isfile(logs_folder + "/scan_" + str(log_file_index) + ".log"):
        log_file_index += 1
    else:
        log_filename = logs_folder + "/scan_" + str(log_file_index) + ".log"
        break
log_file = open(log_filename, "w")

def dual_output(output):
    """ Outputs string to both standard output and the log file. """
    print output,
    log_file.write(output)


### Get start time
start_time = time.time()


### Print header
dual_output("                                                       \n")
dual_output(" ##################################################### \n")
dual_output(" ##                                                 ## \n")
dual_output(" ##  Automated RPV analysis for Snowmass 2013       ## \n")
dual_output(" ##                                                 ## \n")
dual_output(" ##     Written by Mike Saelim - mjs496@cornell.edu ## \n")
dual_output(" ##                                                 ## \n")
dual_output(" ##################################################### \n")
dual_output("                                                       \n")
dual_output("  If directed to proceed with the analysis, the log file will be %s. \n" % log_filename)
dual_output("                                                       \n")
dual_output("  Scan description: \n")
dual_output("     %s \n" % scan_description)
dual_output("                                                       \n")
dual_output(" ------------------------------------------------------------------------- \n")
dual_output("                                                       \n")
dual_output("                                                       \n")




################### Pre-run checking for scan_mode = "background" ####################
if scan_mode == "background":
    # If the sample label already exists, skip the sample
    if os.path.isfile(results_filename):
        for (sample_label, kfactor, lhe_filename) in list(backgrounds_list):
	    results_file = open(results_filename, "r")
            found_sample = False
            for line in results_file:
                if sample_label in line:
                    dual_output("    WARNING: Sample label already present in output file, skipping: %s \n" % sample_label)
                    backgrounds_list.remove((sample_label, kfactor, lhe_filename))
                    break
	    results_file.close()

    dual_output("  There are %d samples to analyze. \n" % len(backgrounds_list))
    dual_output(" \n")


################### Pre-run checking for scan_mode = "signal" ####################
if scan_mode == "signal":
    # If the K-factor file doesn't exist, quit
    if not os.path.isfile(kfactors_filename):
        dual_output("  ERROR: K-factor file %s not found! \n" % kfactors_filename)
        log_file.close()
        os.remove(log_filename)
        exit(0)
        
    # If no K-factor can be found in the K-factor file, skip the point
    for (gluino_mass, stop_mass) in list(signals_list):
        kfactors_file = open(kfactors_filename, "r")
        for line in kfactors_file:
            split_line = line.split()
            if len(split_line) > 0 and split_line[0] == str(gluino_mass):
                break
        else:
            dual_output("    WARNING: No K-factor found for m(gluino) = %d, skipping (%d, %d) \n" % (gluino_mass, gluino_mass, stop_mass))
            signals_list.remove((gluino_mass, stop_mass))
	kfactors_file.close()
        
    # If the sample label already exists, skip the point
    if os.path.isfile(results_filename):
        for (gluino_mass, stop_mass) in list(signals_list):
	    results_file = open(results_filename, "r")
            sample_label = get_sample_label(gluino_mass, stop_mass)
            for line in results_file:
                if sample_label in line:
                    dual_output("    WARNING: Sample label already present in output file, skipping: %s \n" % sample_label)
                    signals_list.remove((gluino_mass, stop_mass))
                    break
	    results_file.close()

    dual_output("  There are %d points to analyze. \n" % len(signals_list))    
    dual_output(" \n")


        
######################################## Confirmation #####################################        

yesno = raw_input("  Would you like to proceed with the analysis? (y/n) ")

if yesno == "y":
    print "  Proceeding..."
    print
else:
    print "  Exiting program..."
    log_file.close()
    os.remove(log_filename)
    exit(0)





######################################## Scanning for scan_mode = "background" #######################################

if scan_mode == "background":
    for (i_loop, (sample_label, kfactor, lhe_filename)) in enumerate(backgrounds_list):
        dual_output("  Loop %d: Running on sample %s, %s ... \n" % (i_loop+1, sample_label, lhe_filename))

        # Unzip the LHE file
        os.system("gunzip " + lhe_filename + " &> /dev/null")
        
        # Run the Hadronizer
        # Catch fails in case Pythia aborts or segfaults
        for i_attempt in range(max_attempts_per_sample):
            if os.path.isfile(hepmc_file_folder + "/" + sample_label + ".hepmc"):
                os.remove(hepmc_file_folder + "/" + sample_label + ".hepmc")
            if os.system(hadronizer_exec + " " + hepmc_file_folder + "/" + sample_label + ".hepmc LHE " + lhe_filename[:-3] + " &> /dev/null") == 0:
                dual_output("    Hadronizer done! \n")
                break
            dual_output("    Hadronizer failed... \n")
        else:
            dual_output("    ERROR: Hadronizer failed %d times, skipping this sample! \n" % max_attempts_per_sample)
            error_count += 1
            continue

        # Rezip the LHE file
        os.system("gzip " + lhe_filename + " &> /dev/null")

        # Run Delphes
        # Catch fails in case Delphes aborts or segfaults
        for i_attempt in range(max_attempts_per_sample):
            if os.path.isfile(root_file_folder + "/" + sample_label + ".root"):
                os.remove(root_file_folder + "/" + sample_label + ".root")
            if os.system(delphes_exec + " " + delphes_card_file + " " + \
                         root_file_folder + "/" + sample_label + ".root " + \
                         hepmc_file_folder + "/" + sample_label + ".hepmc &> /dev/null") == 0:
                dual_output("    Delphes done! \n")
                break
            dual_output("    Delphes failed... \n")
        else:
            dual_output("    ERROR: Delphes failed %d times, skipping this sample! \n" % max_attempts_per_sample)
            error_count += 1
            continue
        
        # Extract LO cross section (in pb) from HepMC file
        # The last entry in the file should have the most accurate cross section
        hepmc_file = open(hepmc_file_folder + "/" + sample_label + ".hepmc", "r")
        lo_cross_section = -1.
        for line in hepmc_file:
            split_line = line.split()
            if len(split_line) > 0 and split_line[0] == "C":
                lo_cross_section = split_line[1]
        hepmc_file.close()

	# Remove HepMC file
	os.remove(hepmc_file_folder + "/" + sample_label + ".hepmc")

        # Run analysis
        os.system(analysis_exec + " " + results_filename + " '" + sample_label + "' " + str(lo_cross_section) + " " + \
                  str(kfactor) + " " + root_file_folder + "/" + sample_label + ".root &> /dev/null")
        dual_output("    Analysis done! \n")




######################################## Scanning for scan_mode = "signal" #######################################

if scan_mode == "signal":
    for (i_loop, (gluino_mass, stop_mass)) in enumerate(signals_list):
        sample_label = get_sample_label(gluino_mass, stop_mass)
        dual_output("  Loop %d: Running on sample %s ... \n" % (i_loop+1, sample_label))

        # Generate input .spc and .cmnd files for this point
        os.system(geninputs_exec + " rpv_gluino_stop.spc rpv_gluinos.cmnd " + str(num_events) + " 14 " + \
                  str(gluino_mass) + " " + str(stop_mass) + " &> /dev/null")
        dual_output("    Inputs generated! \n")
        
        # Run the Hadronizer
        # Catch fails in case Pythia aborts or segfaults
        for i_attempt in range(max_attempts_per_sample):
            if os.path.isfile(hepmc_file_folder + "/" + sample_label + ".hepmc"):
                os.remove(hepmc_file_folder + "/" + sample_label + ".hepmc")
            if os.system(hadronizer_exec + " " + hepmc_file_folder + "/" + sample_label + ".hepmc CMND rpv_gluinos.cmnd rpv_gluino_stop.spc &> /dev/null") == 0:
                dual_output("    Hadronizer done! \n")
                break
            dual_output("    Hadronizer failed... \n")
        else:
            dual_output("    ERROR: Hadronizer failed %d times, skipping this sample! \n" % max_attempts_per_sample)
            error_count += 1
            continue

        # Run Delphes
        # Catch fails in case Delphes aborts or segfaults
        for i_attempt in range(max_attempts_per_sample):
            if os.path.isfile(root_file_folder + "/" + sample_label + ".root"):
                os.remove(root_file_folder + "/" + sample_label + ".root")
            if os.system(delphes_exec + " " + delphes_card_file + " " + \
                         root_file_folder + "/" + sample_label + ".root " + \
                         hepmc_file_folder + "/" + sample_label + ".hepmc &> /dev/null") == 0:
                dual_output("    Delphes done! \n")
                break
            dual_output("    Delphes failed... \n")
        else:
            dual_output("    ERROR: Delphes failed %d times, skipping this sample! \n" % max_attempts_per_sample)
            error_count += 1
            continue
        
        # Extract LO cross section (in pb) from HepMC file
        # The last entry in the file should have the most accurate cross section
        hepmc_file = open(hepmc_file_folder + "/" + sample_label + ".hepmc", "r")
        lo_cross_section = -1.
        for line in hepmc_file:
            split_line = line.split()
            if len(split_line) > 0 and split_line[0] == "C":
                lo_cross_section = split_line[1]
        hepmc_file.close()

	# Remove HepMC file
	os.remove(hepmc_file_folder + "/" + sample_label + ".hepmc")

        # Extract K-factor from K-factor file
        kfactors_file = open(kfactors_filename, "r")
        kfactor = -1.
        for line in kfactors_file:
            split_line = line.split()
            if len(split_line) > 0 and split_line[0] == str(gluino_mass):
                kfactor = split_line[1]
                break
        kfactors_file.close()

        # Run analysis
        os.system(analysis_exec + " " + results_filename + " '" + sample_label + "' " + str(lo_cross_section) + " " + \
                  str(kfactor) + " " + root_file_folder + "/" + sample_label + ".root &> /dev/null")
        dual_output("    Analysis done! \n")




########################################## Clean up and output final statistics ####################################

end_time = time.time()
elapsed_time = end_time - start_time

dual_output("                                                        \n")
dual_output("                                                        \n")
dual_output("  Scan finished!                                        \n")
if scan_mode == "background":
    dual_output("    Analyzed %d samples                             \n" % len(backgrounds_list))
if scan_mode == "signal":
    dual_output("    Analyzed %d samples                             \n" % len(signals_list))
dual_output("    Elapsed wall clock time:  %d days, %02d:%02d:%02d   \n" % (int(elapsed_time) / (24*60*60), \
                                                                            int(elapsed_time) % (24*60*60) / (60*60), \
                                                                            int(elapsed_time) % (60*60) / 60, \
                                                                            int(elapsed_time) % 60))
dual_output("                                                        \n")

if error_count == 0:
    dual_output("  There were no errors!  Hooray!                    \n")
else:
    dual_output("  There were %d errors - go check them, dummy!      \n" % error_count)

dual_output("                                                        \n")



log_file.close()















