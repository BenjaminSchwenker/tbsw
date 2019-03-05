// SeedGenerator implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// Local includes
#include "SeedGenerator.h"
#include "StraightLineTrackModel.h"

// C++ includes
#include <iostream>
#include <cassert>
#include <cstdlib>



using namespace std;

namespace depfet {

SeedGenerator::SeedGenerator( double acharge, double amom ) : charge(acharge), mom(amom) {;}

// Create a seed track state 
TBTrackState SeedGenerator::CreateSeedTrack(const TBHit& FirstHit, const TBHit& SecondHit, const TBDetector& Detector)  const
{  

  // Check hits are on different detectors 
  if ( FirstHit.GetSensorID() == SecondHit.GetSensorID() ) {
    return CreateSeedTrack(FirstHit, Detector);  
  }     

  // Get plane numbers 
  int firstplane = Detector.GetPlaneNumber( FirstHit.GetSensorID() );
  int secondplane = Detector.GetPlaneNumber( SecondHit.GetSensorID() );
  
  // Sort hits along beam line 
  if ( firstplane > secondplane ) {
    // just call method with swapped order of hits
    return CreateSeedTrack(SecondHit, FirstHit, Detector); 
  } 
  
  // Compute global space points
  const ReferenceFrame& FirstFrame = Detector.GetDet( firstplane ).GetNominal(); 
  const ReferenceFrame& SecondFrame = Detector.GetDet( secondplane ).GetNominal(); 
   
  Vector3d FirstPoint = FirstHit.GetLocalSpacePoint(); 
  Vector3d SecondPoint = SecondHit.GetLocalSpacePoint();
   
  // Compute local track direction 
  Vector3d GobalDirection = SecondFrame.TransformPointToGlobal( SecondPoint ) - FirstFrame.TransformPointToGlobal( FirstPoint );
  Vector3d LocalDirection = FirstFrame.TransformVecToLocal( GobalDirection );  
  
  // Seed parameters at first sensor
  TrackState Pars;  
  Pars[0] = LocalDirection[0]/LocalDirection[2]; 
  Pars[1] = LocalDirection[1]/LocalDirection[2];
  Pars[2] = FirstHit.GetCoord()[0];
  Pars[3] = FirstHit.GetCoord()[1];
  Pars[4] = charge/mom;
  
  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane(firstplane); 
   
  return Seed; 
}

// Create a seed track state  
TBTrackState SeedGenerator::CreateSeedTrack(const TBHit& Hit, const TBDetector& Detector) const
{ 

  // Get plane number of hit
  int planenumber = Detector.GetPlaneNumber( Hit.GetSensorID() );  

  // Compute global space point
  const ReferenceFrame& Frame = Detector.GetDet( planenumber ).GetNominal(); 
  
  // Seed follows z direction
  Vector3d GobalDirection;
  GobalDirection << 0,0,1;   

  // Compute local track direction 
  Vector3d LocalDirection = Frame.TransformVecToLocal(GobalDirection);  
  
  TrackState Pars; 
  Pars[0] = LocalDirection[0]/LocalDirection[2];  
  Pars[1] = LocalDirection[1]/LocalDirection[2];
  Pars[2] = Hit.GetCoord()[0];
  Pars[3] = Hit.GetCoord()[1];
  Pars[4] = charge/mom;

  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane(planenumber); 
        
  return Seed; 
}


// Create a seed track state  
TBTrackState SeedGenerator::CreateSeedTrack() const
{ 
 
  // Seed crosses origin and follows z axis 
  TrackState Pars; 
  Pars[0] = 0; 
  Pars[1] = 0; 
  Pars[2] = 0;
  Pars[3] = 0;
  Pars[4] = charge/mom;
        
  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane(0); 
        
  return Seed; 
}



 
} // Namespace;
