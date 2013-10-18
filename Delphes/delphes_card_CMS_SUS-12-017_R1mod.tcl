#######################################################################################
###   delphes_card_CMS-SUS-12-017_R1mod.tcl
###
###   Written by: Mike Saelim, mjs496@cornell.edu, 9/3/13
###
###
###
### This Delphes card was written for the Snowmass 2013 implementation of our analysis in arXiv:1302.2146.
###  - It handles detector simulation and object reconstruction, but leaves the event-level cuts and analysis for later code.
###
### We implement the CMS SSDL+b jets+MET search at 8 TeV and 10.5 fb^-1 (CMS PAS-SUS-12-029/CMS SUS-12-017), but with jet R = 1.0.
###  - The leptonic analysis is our special isolation-aware modification of their suggested analysis.
###  - We cluster jets with R = 1.0 anti-kt.
###  - The b-tagging is done with jet pT-dependent tagging and mistagging efficiencies for the CSVM tagger, from arXiv:1211.4462/BTV-12-001.
###  - We ignore the H_T and MET turn-on efficiencies provided by the paper and define them properly.
###
### It is based off delphes_card_CMS.tcl included in the Delphes /examples/ folder.
### Major changes:
###  - removed photon reconstruction, GenJet reconstruction, and tau tagging entirely
###  - modified Efficiency module (SaelimEfficiency) to handle ElectronEfficiency and MuonEfficiency, allowing for two working points
###       (this way, leptons passing loose ID requirements automatically pass tight ID requirements)
###  - separate Isolation modules for loose and tight lepton selections
###  - included lepton ID pT and eta cuts, and ID efficiencies, in the *Efficiency modules
###  - included lepton isolation cuts in the *IsolationTight and *IsolationLoose modules
###  - included jet R and pT cut in FastJetFinder module
###  - included jet eta cut in JetEnergyScale module
###  - included jet R, b-tagging, and mistag efficiencies in BTagging module
###  - included only jets in the ScalarHT module
#######################################################################################




#######################################
# Order of execution of various modules
#######################################
set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing
  
  TrackMerger
  Calorimeter
  EFlowMerger

  ElectronEfficiency
  ElectronIsolationTight
  ElectronIsolationLoose

  MuonEfficiency
  MuonIsolationTight
  MuonIsolationLoose

  FastJetFinder
  JetEnergyScale
  BTagging

  UniqueObjectFinder
  ScalarHT
  MissingET

  TreeWriter
}



#################################
# Propagate particles in cylinder
#################################
module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.00

  # magnetic field
  set Bz 3.8
}



####################################
# Charged hadron tracking efficiency
####################################
module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.60) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.85) + \
                         (abs(eta) > 2.5)                                                  * (0.00)}
}



##############################
# Electron tracking efficiency
##############################
module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.73) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2) * (0.95) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.99) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.50) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.83) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.90) + \
                         (abs(eta) > 2.5)                                                  * (0.00)}
}



##########################
# Muon tracking efficiency
##########################
module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.99) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.98) + \
                         (abs(eta) > 2.5)                                                  * (0.00)}
}



########################################
# Momentum resolution for charged tracks
########################################
module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0e1) * (0.20) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.20) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.20) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.20) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.20) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.20)}
}



#################################
# Energy resolution for electrons
#################################
module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons
  set ResolutionFormula {                  (abs(eta) <= 2.5) * (energy > 0.1   && energy <= 2.0e1) * (0.0225)*energy + \
                                           (abs(eta) <= 2.5) * (energy > 2.0e1)                    * sqrt(energy^2*0.007^2 + energy*0.07^2 + 0.35^2) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0)                                       * sqrt(energy^2*0.007^2 + energy*0.07^2 + 0.35^2) + \
                         (abs(eta) > 3.0 && abs(eta) <= 5.0)                                       * sqrt(energy^2*0.107^2 + energy*2.08^2)}

}



###############################
# Momentum resolution for muons
###############################
module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for muons
  set ResolutionFormula {                  (abs(eta) <= 0.5) * (pt > 0.1   && pt <= 5.0)   * (0.02) + \
                                           (abs(eta) <= 0.5) * (pt > 5.0   && pt <= 1.0e2) * (0.015) + \
                                           (abs(eta) <= 0.5) * (pt > 1.0e2)                * (0.03) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1   && pt <= 5.0)   * (0.03) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 5.0   && pt <= 1.0e2) * (0.02) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 5.0)   * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 5.0   && pt <= 1.0e2) * (0.035) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.05)}
}



##############
# Track merger
##############
module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  set OutputArray tracks
}



#############
# Calorimeter
#############
module Calorimeter Calorimeter {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray towers
  set PhotonOutputArray photons

  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowTowers

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 5 degrees towers
  set PhiBins {}
  for {set i -36} {$i <= 36} {incr i} {
    add PhiBins [expr {$i * $pi/36.0}]
  }
  foreach eta {-1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.653} {
    add EtaPhiBins $eta $PhiBins
  }

  # 10 degrees towers
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }
  foreach eta {-4.35 -4.175 -4 -3.825 -3.65 -3.475 -3.3 -3.125 -2.95 -2.868 -2.65 -2.5 -2.322 -2.172 -2.043 -1.93 -1.83 -1.74 -1.653 1.74 1.83 1.93 2.043 2.172 2.322 2.5 2.65 2.868 2.95 3.125 3.3 3.475 3.65 3.825 4 4.175 4.35 4.525} {
    add EtaPhiBins $eta $PhiBins
  }

  # 20 degrees towers
  set PhiBins {}
  for {set i -9} {$i <= 9} {incr i} {
    add PhiBins [expr {$i * $pi/9.0}]
  }
  foreach eta {-5 -4.7 -4.525 4.7 5} {
    add EtaPhiBins $eta $PhiBins
  }

  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {0.0 1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0 0.0}
  add EnergyFraction {22} {1.0 0.0}
  add EnergyFraction {111} {1.0 0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0 0.0}
  add EnergyFraction {13} {0.0 0.0}
  add EnergyFraction {14} {0.0 0.0}
  add EnergyFraction {16} {0.0 0.0}
  add EnergyFraction {1000022} {0.0 0.0}
  add EnergyFraction {1000023} {0.0 0.0}
  add EnergyFraction {1000025} {0.0 0.0}
  add EnergyFraction {1000035} {0.0 0.0}
  add EnergyFraction {1000045} {0.0 0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3 0.7}
  add EnergyFraction {3122} {0.3 0.7}

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}
  set ECalResolutionFormula {                  (abs(eta) <= 3.0) * sqrt(energy^2*0.007^2 + energy*0.07^2 + 0.35^2)  + \
                             (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.107^2 + energy*2.08^2)}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set HCalResolutionFormula {                  (abs(eta) <= 1.5) * sqrt(energy^2*0.055^2 + energy*1.20^2 + 3.50^2) + \
                                               (abs(eta) > 1.5 && abs(eta) <= 3)   * sqrt(energy^2*0.065^2 + energy*2.00^2 + 3.50^2) + \
                                               (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.11^2 + energy*3.13^2)}
}



####################
# Energy flow merger
####################
module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray Calorimeter/eflowTracks
  add InputArray Calorimeter/eflowTowers
  add InputArray MuonMomentumSmearing/muons
  set OutputArray eflow
}



#####################
# Electron efficiency
#####################
module SaelimEfficiency ElectronEfficiency {
  set InputArray ElectronEnergySmearing/electrons
  set OutputArray1 tightElectrons
  set OutputArray2 looseElectrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
  set EfficiencyFormula1 {                  (abs(eta) <= 1.442) * (pt > 20.0) * (0.73) + \
			   (abs(eta) > 1.566 && abs(eta) < 2.4) * (pt > 20.0) * (0.73)}
  set EfficiencyFormula2 {                  (abs(eta) <= 1.442) * (pt > 5.0)  * (0.73) + \
			   (abs(eta) > 1.566 && abs(eta) < 2.4) * (pt > 5.0)  * (0.73)}
}



####################
# Electron isolation
####################
module Isolation ElectronIsolationTight {
  set CandidateInputArray ElectronEfficiency/tightElectrons
  set IsolationInputArray EFlowMerger/eflow

  set OutputArray tightElectrons

  set DeltaRMax 0.3

  set PTMin 0.0

  set PTRatioMax 0.1
}
module Isolation ElectronIsolationLoose {
  set CandidateInputArray ElectronEfficiency/looseElectrons
  set IsolationInputArray EFlowMerger/eflow

  set OutputArray looseElectrons

  set DeltaRMax 0.3

  set PTMin 0.0

  set PTRatioMax 0.2
}



#################
# Muon efficiency
#################
module SaelimEfficiency MuonEfficiency {
  set InputArray MuonMomentumSmearing/muons
  set OutputArray1 tightMuons
  set OutputArray2 looseMuons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
  set EfficiencyFormula1 { (abs(eta) <= 2.4) * (pt > 20.0) * (0.84)}
  set EfficiencyFormula2 { (abs(eta) <= 2.4) * (pt > 5.0)  * (0.84)}
}



################
# Muon isolation
################
module Isolation MuonIsolationTight {
  set CandidateInputArray MuonEfficiency/tightMuons
  set IsolationInputArray EFlowMerger/eflow

  set OutputArray tightMuons

  set DeltaRMax 0.3

  set PTMin 0.0

  set PTRatioMax 0.1
}
module Isolation MuonIsolationLoose {
  set CandidateInputArray MuonEfficiency/looseMuons
  set IsolationInputArray EFlowMerger/eflow

  set OutputArray looseMuons

  set DeltaRMax 0.3

  set PTMin 0.0

  set PTRatioMax 0.2
}



############
# Jet finder
############
module FastJetFinder FastJetFinder {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.0

  set JetPTMin 40.0
}



##################
# Jet Energy Scale
##################
module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula { (abs(eta) <= 2.4) * (1.08) + \
		     (abs(eta) >  2.4) * (0.0)}
}



###########
# b-tagging
###########
module BTagging BTagging {
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  set DeltaR 1.0

  # This seems to need to be set to a nonzero value, or else you get a warning:
  # "Warning in <TVector3::PseudoRapidity>: transvers momentum = 0! return +/- 10e10"
  set PartonPTMin 1.0

  set PartonEtaMax 5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # BTV-12-001
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {               (pt <=  10.4) * (0.0132) + (pt >  10.4 && pt <=  20.7) * (0.0176) + \ 
                              (pt >  20.7 && pt <=  31.1) * (0.0119) + (pt >  31.1 && pt <=  41.4) * (0.0147) + \ 
                              (pt >  41.4 && pt <=  51.8) * (0.0146) + (pt >  51.8 && pt <=  62.1) * (0.0135) + \
                              (pt >  62.1 && pt <=  72.5) * (0.0153) + (pt >  72.5 && pt <=  82.8) * (0.0136) + \
                              (pt >  82.8 && pt <=  93.2) * (0.0159) + (pt >  93.2 && pt <= 103.6) * (0.0154) + \
                              (pt > 103.6 && pt <= 113.9) * (0.0170) + (pt > 113.9 && pt <= 124.3) * (0.0174) + \
                              (pt > 124.3 && pt <= 134.6) * (0.0177) + (pt > 134.6 && pt <= 145.0) * (0.0178) + \
                              (pt > 145.0 && pt <= 155.3) * (0.0199) + (pt > 155.3 && pt <= 166.0) * (0.0191) + \
                              (pt > 166.0 && pt <= 176.4) * (0.0198) + (pt > 176.4 && pt <= 186.4) * (0.0202) + \
                              (pt > 186.4 && pt <= 196.8) * (0.0213) + (pt > 196.8 && pt <= 207.1) * (0.0206) + \ 
                              (pt > 207.1 && pt <= 217.5) * (0.0206) + (pt > 217.5 && pt <= 228.2) * (0.0222) + \ 
                              (pt > 228.2 && pt <= 238.5) * (0.0233) + (pt > 238.5 && pt <= 248.5) * (0.0233) + \ 
                              (pt > 248.5 && pt <= 258.9) * (0.0221) + (pt > 258.9 && pt <= 269.3) * (0.0239) + \ 
                              (pt > 269.3 && pt <= 279.6) * (0.0244) + (pt > 279.6 && pt <= 290.0) * (0.0244) + \ 
                              (pt > 290.0 && pt <= 300.3) * (0.0246) + (pt > 300.3 && pt <= 310.7) * (0.0260) + \ 
                              (pt > 310.7 && pt <= 321.0) * (0.0243) + (pt > 321.0 && pt <= 331.4) * (0.0253) + \ 
                              (pt > 331.4 && pt <= 341.7) * (0.0254) + (pt > 341.7 && pt <= 352.1) * (0.0268) + \ 
                              (pt > 352.1 && pt <= 362.5) * (0.0257) + (pt > 362.5 && pt <= 372.8) * (0.0270) + \ 
                              (pt > 372.8 && pt <= 383.2) * (0.0274) + (pt > 383.2 && pt <= 393.5) * (0.0270) + \ 
                              (pt > 393.5 && pt <= 403.9) * (0.0288) + (pt > 403.9 && pt <= 414.2) * (0.0299) + \ 
                              (pt > 414.2 && pt <= 424.6) * (0.0292) + (pt > 424.6 && pt <= 435.0) * (0.0289) + \ 
                              (pt > 435.0 && pt <= 445.3) * (0.0302) + (pt > 445.3 && pt <= 455.7) * (0.0298) + \ 
                              (pt > 455.7 && pt <= 466.0) * (0.0310) + (pt > 466.0 && pt <= 476.4) * (0.0305) + \ 
                              (pt > 476.4 && pt <= 486.7) * (0.0308) + (pt > 486.7 && pt <= 497.1) * (0.0299) + \ 
                              (pt > 497.1 && pt <= 507.4) * (0.0315) + (pt > 507.4 && pt <= 517.8) * (0.0319) + \ 
                              (pt > 517.8 && pt <= 528.2) * (0.0344) + (pt > 528.2 && pt <= 538.5) * (0.0324) + \ 
                              (pt > 538.5 && pt <= 548.9) * (0.0340) + (pt > 548.9 && pt <= 559.2) * (0.0315) + \ 
                              (pt > 559.2 && pt <= 569.6) * (0.0322) + (pt > 569.6 && pt <= 579.9) * (0.0333) + \ 
                              (pt > 579.9 && pt <= 590.3) * (0.0301) + (pt > 590.3 && pt <= 600.6) * (0.0317) + \ 
                              (pt > 600.6 && pt <= 611.0) * (0.0317) + (pt > 611.0 && pt <= 621.4) * (0.0324) + \ 
                              (pt > 621.4 && pt <= 631.7) * (0.0344) + (pt > 631.7 && pt <= 642.1) * (0.0364) + \ 
                              (pt > 642.1 && pt <= 652.4) * (0.0315) + (pt > 652.4 && pt <= 662.5) * (0.0329) + \ 
                              (pt > 662.5)                * (0.0334)} 

  # efficiency formula for b-jets
  add EfficiencyFormula {5} {               (pt <=  40.0) * (0.000) + \
                              (pt >  40.0 && pt <=  50.0) * (0.640) + \
                              (pt >  50.0 && pt <=  60.0) * (0.704) + \
                              (pt >  60.0 && pt <=  70.0) * (0.710) + \
                              (pt >  70.0 && pt <=  80.0) * (0.725) + \
                              (pt >  80.0 && pt <= 100.0) * (0.750) + \
                              (pt > 100.0 && pt <= 120.0) * (0.752) + \
                              (pt > 120.0 && pt <= 160.0) * (0.752) + \
                              (pt > 160.0 && pt <= 210.0) * (0.732) + \
                              (pt > 210.0 && pt <= 260.0) * (0.712) + \
                              (pt > 260.0 && pt <= 320.0) * (0.707) + \
                              (pt > 320.0 && pt <= 400.0) * (0.687) + \
                              (pt > 400.0 && pt <= 500.0) * (0.668) + \
                              (pt > 500.0)                * (0.580)}
}



#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################
module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray ElectronIsolationTight/tightElectrons electrons
  add InputArray MuonIsolationTight/tightMuons muons
  add InputArray JetEnergyScale/jets jets
}



##################
# Scalar HT merger
##################
module Merger ScalarHT {
# add InputArray InputArray
  add InputArray UniqueObjectFinder/jets
  set EnergyOutputArray energy
}



###################
# Missing ET merger
###################
module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}



##################
# ROOT tree writer
##################
module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/muons Muon Muon
  add Branch ScalarHT/energy ScalarHT ScalarHT
  add Branch MissingET/momentum MissingET MissingET
  add Branch ElectronIsolationLoose/looseElectrons looseElectron Electron
  add Branch MuonIsolationLoose/looseMuons looseMuon Muon
}

