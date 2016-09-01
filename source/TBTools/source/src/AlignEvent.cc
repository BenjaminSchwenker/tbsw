#include "AlignEvent.h"
 	 	
using namespace std;

ClassImp(depfet::AlignEvent)
 	
namespace depfet {	

AlignEvent::AlignEvent()
{
  fIndex = new TArrayI;
  fMeasurements = new TVectorD;
  fMeasuredCovariance = new TMatrixDSym;
  fRefTrackParameters = new TVectorD;
}
 	
AlignEvent::~AlignEvent()
{
  delete fIndex;
  delete fMeasurements;
  delete fMeasuredCovariance;
  delete fRefTrackParameters; 
}

} // Namespace;
