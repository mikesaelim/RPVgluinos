
/** \class SaelimEfficiency
 *
 *  Selects candidates from the InputArray according to two efficiency formulas.
 *
 *  modified from Efficiency.cc on 9/3/13 by Mike Saelim, mjs496@cornell.edu
 *
 *
 *  $Date: 2013-02-12 14:57:44 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 905 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/SaelimEfficiency.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

SaelimEfficiency::SaelimEfficiency() :
  fFormula1(0), fFormula2(0), fItInputArray(0)
{
  fFormula1 = new DelphesFormula;
  fFormula2 = new DelphesFormula;
}

//------------------------------------------------------------------------------

SaelimEfficiency::~SaelimEfficiency()
{
  if(fFormula1) delete fFormula1;
  if(fFormula2) delete fFormula2;
}

//------------------------------------------------------------------------------

void SaelimEfficiency::Init()
{
  // read efficiency formula

  fFormula1->Compile(GetString("EfficiencyFormula1", "1.0"));
  fFormula2->Compile(GetString("EfficiencyFormula2", "1.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray1 = ExportArray(GetString("OutputArray1", "stableParticles1"));
  fOutputArray2 = ExportArray(GetString("OutputArray2", "stableParticles2"));
}

//------------------------------------------------------------------------------

void SaelimEfficiency::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void SaelimEfficiency::Process()
{ 
  Candidate *candidate;
  Double_t pt, eta, phi;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();

    // apply efficiency formulas
    Double_t random_number = gRandom->Uniform();
    if(random_number <= fFormula1->Eval(pt, eta))  fOutputArray1->Add(candidate);
    if(random_number <= fFormula2->Eval(pt, eta))  fOutputArray2->Add(candidate);
  }
}

//------------------------------------------------------------------------------
