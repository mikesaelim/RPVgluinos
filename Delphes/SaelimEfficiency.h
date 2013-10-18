#ifndef SaelimEfficiency_h
#define SaelimEfficiency_h

/** \class SaelimEfficiency
 *
 *  Selects candidates from the InputArray according to two efficiency formulas.
 *
 *  modified from Efficiency.h on 9/3/13 by Mike Saelim, mjs496@cornell.edu
 *
 *
 *  $Date: 2013-02-12 14:57:44 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 905 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class SaelimEfficiency: public DelphesModule
{
public:

  SaelimEfficiency();
  ~SaelimEfficiency();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula1; //!
  DelphesFormula *fFormula2;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray1; //!
  TObjArray *fOutputArray2;

  ClassDef(SaelimEfficiency, 1)
};

#endif
