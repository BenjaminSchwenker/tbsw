// TBTrackElement implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// TBTools includes 
#include "TBTrackElement.h"

using namespace std;

namespace depfet {


/* Constructor 
 */
TBTrackElement::TBTrackElement(Det& aDetUnit) : DetUnit(aDetUnit) 
{
  CrossedFlag = false; 
  LocalChiSqu = 0;
  hasHit=false;
}

/** Set measured hit
 */
void TBTrackElement::SetHit(const TBHit& aHit)
{
  if(hasHit){
      HitStore[0]=aHit;
  } else{
    HitStore.push_back(aHit);
  }
  hasHit=true;

}
  
/** Remove (bad) hit - Outlier rejection  
 */
void TBTrackElement::RemoveHit() {
  HitStore.clear()
  hasHit=false;
}

/** Get selected hit - query HasHit() befor 
 */
TBHit& TBTrackElement::GetHit()
{
   if(!hasHit){
       throw std::out_of_range("Track Element has no hit!");
   }
   return HitStore[0];
}

/** Set track state
   */
void TBTrackElement::SetState(const TBTrackState& aState)
{
  State = aState;
}

} // Namespace;


