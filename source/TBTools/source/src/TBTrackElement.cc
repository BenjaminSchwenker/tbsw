// TBTrackElement implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// DEPFETTrackTools includes 
#include "TBTrackElement.h"

using namespace CLHEP;
using namespace std;

namespace depfet {


/* Constructor 
 */
TBTrackElement::TBTrackElement(Det& aDetUnit) : DetUnit(aDetUnit) 
{
  CrossedFlag = false; 
  LocalChiSqu = 0;   
}

/** Set measured hit
 */
void TBTrackElement::SetHit(TBHit& aHit)
{
  if (HitStore.size() == 0) HitStore.push_back(aHit); 
  else HitStore[0] = aHit; 
}
  
/** True if hit is set  
 */
bool TBTrackElement::HasHit()
{
  return (HitStore.size() == 1); 
} 
  
/** Remove (bad) hit - Outlier rejection  
 */
void TBTrackElement::RemoveHit() {
  HitStore.clear();  
}

/** Get selected hit - query HasHit() befor 
 */
TBHit& TBTrackElement::GetHit()
{
  return HitStore[0];
}

/** Set track state
   */
void TBTrackElement::SetState(TBTrackState& aState) 
{
  State = aState;
}

} // Namespace;


