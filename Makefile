CPP = g++
CXXFLAGS = -O2 -ansi -pedantic -Wall -Wextra -Wshadow -fbounds-check
PYTHIA_LIB = /Users/donerkebab/Desktop/work/software/pythia8176/lib/archive
PYTHIA_INC = /Users/donerkebab/Desktop/work/software/pythia8176/include
FASTJETLOCATION = /usr/local
HEPMCLOCATION = /Users/donerkebab/Desktop/work/software/HepMC/x86_64-mac106-gcc42-opt
DELPHESLOCATION = /Users/donerkebab/Desktop/work/software/Delphes-3.0.10
ROOTLOCATION = /Users/donerkebab/Desktop/work/software/root


GenerateInputs: 
	$(CPP) -I./ $@.cc $(CXXFLAGS) -Wno-shadow -o $@

Hadronizer: $(PYTHIA_LIB)/libpythia8.a $(PYTHIA_LIB)/libhepmcinterface.a
	$(CPP) -I$(PYTHIA_INC) -I$(HEPMCLOCATION)/include -I./ $@.cc $(CXXFLAGS) -Wno-shadow -o $@ -L$(PYTHIA_LIB) -L$(HEPMCLOCATION)/lib -lpythia8 -llhapdfdummy -lhepmcinterface -lhepmc

JetMassCutAnalysis:
	$(CPP) -I./ -I$(DELPHESLOCATION) -I$(ROOTLOCATION)/include $@.cc `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` $(CXXFLAGS) -Wno-shadow -Wno-long-long -o $@ `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` `$(ROOTLOCATION)/bin/root-config --cflags --glibs` -L$(DELPHESLOCATION) -lDelphes




.PHONY: GenerateInputs Hadronizer JetMassCutAnalysis