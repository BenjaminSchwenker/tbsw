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

// CLHEP includes 
#include <CLHEP/Matrix/Matrix.h>

using namespace std;
using namespace CLHEP;

namespace depfet {

SeedGenerator::SeedGenerator( double acharge, double amom ) : charge(acharge), mom(amom) {;}

// Create a seed track state 
TBTrackState SeedGenerator::CreateSeedTrack(TBHit FirstHit, TBHit SecondHit, TBDetector& Detector)   
{  

  // Check hits are on different detectors 
  if ( FirstHit.GetDAQID() == SecondHit.GetDAQID() ) {
    return CreateSeedTrack(FirstHit, Detector);  
  }     

  // Get plane numbers 
  int firstplane = Detector.GetPlaneNumber( FirstHit.GetDAQID() );
  int secondplane = Detector.GetPlaneNumber( SecondHit.GetDAQID() );
  
  // Sort hits along beam line 
  if ( firstplane > secondplane ) {
    TBHit TmpHit = FirstHit; 
    FirstHit = SecondHit; 
    SecondHit = TmpHit; 
  } 
  
  // Compute global space points
  ReferenceFrame FirstFrame = Detector.GetDet( firstplane ).GetNominal(); 
  ReferenceFrame SecondFrame = Detector.GetDet( secondplane ).GetNominal(); 
   
  HepVector FirstPoint = FirstHit.GetLocalSpacePoint(); 
  HepVector FirstGPoint = FirstFrame.TransformPointToGlobal(FirstPoint);
  
  HepVector SecondPoint = SecondHit.GetLocalSpacePoint(); 
  HepVector SecondGPoint = SecondFrame.TransformPointToGlobal(SecondPoint);
  
  // Compute global track direction 
  HepVector GobalDirection = SecondGPoint - FirstGPoint;
 
  // Compute local track direction 
  HepVector LocalDirection = FirstFrame.TransformVecToLocal(GobalDirection);  
  
  // Seed parameters at first sensor
  HepMatrix Pars(5,1,0);  
  Pars[0][0] = LocalDirection[0]/LocalDirection[2]; 
  Pars[1][0] = LocalDirection[1]/LocalDirection[2];
  Pars[2][0] = FirstHit.GetCoord()[0][0];
  Pars[3][0] = FirstHit.GetCoord()[1][0];
  Pars[4][0] = charge/mom;
  
  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane(firstplane); 
   
  return Seed; 
}

// Create a seed track state  
TBTrackState SeedGenerator::CreateSeedTrack(TBHit Hit, TBDetector& Detector)
{ 

  // Get plane number of hit
  int planenumber = Detector.GetPlaneNumber( Hit.GetDAQID() );  

  // Compute global space point
  ReferenceFrame Frame = Detector.GetDet( planenumber ).GetNominal(); 
  
  // Seed follows z direction
  HepVector GobalDirection(3,0);
  GobalDirection[2] = 1; 

  // Compute local track direction 
  HepVector LocalDirection = Frame.TransformVecToLocal(GobalDirection);  
  
  HepMatrix Pars(5,1,0); 
  Pars[0][0] = LocalDirection[0]/LocalDirection[2];  
  Pars[1][0] = LocalDirection[1]/LocalDirection[2];
  Pars[2][0] = Hit.GetCoord()[0][0];
  Pars[3][0] = Hit.GetCoord()[1][0];
  Pars[4][0] = charge/mom;

  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane(planenumber); 
        
  return Seed; 
}


// Create a seed track state  
TBTrackState SeedGenerator::CreateSeedTrack()
{ 
 
  
  
  // Seed crosses origin and follows z axis 
  HepMatrix Pars(5,1,0); 
  Pars[0][0] = 0; 
  Pars[1][0] = 0; 
  Pars[2][0] = 0;
  Pars[3][0] = 0;
  Pars[4][0] = charge/mom;
        
  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane(0); 
        
  return Seed; 
}



 
} // Namespace;
