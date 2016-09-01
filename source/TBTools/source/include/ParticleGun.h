#ifndef ParticleGun_H
#define ParticleGun_H 1

// Include DEPFETTrackTools header files
#include "TBDetector.h"
#include "TBHit.h"
#include "TBTrack.h"

#include <TRandom.h>
#include <TRandom3.h>

namespace depfet { 


/** Class ParticleGun
 *  
 *  The ParticleGun can be used for toy simulations of a 
 *  pixel tracking telescope. The ParticleGun emulates 
 *  the properties of the particle beam, multiple scattering 
 *  of tracks in the telescope layers and air and digitization
 *  in pixel sensors. 
 *  
 *  The ParticleGun class provides a user interface to 
 *  specify the relevant properties of the particle beam 
 *  to be simulated. The class has member functions to 
 *  simulate individual particle trajectories and to 
 *  simulate the detector response, i.e. pixel hits.  
 *  
 *  In short: All particles enter the tracking area from a final 
 *  collimator opening of the beam line. The tracker volume is 
 *  filled with air and free of a magnetic field. The tracking 
 *  telescope is specified by a TBDector object. Sensor positions, 
 *  material budget per sensor and spatial resolutions are read 
 *  from the TBDetector object. The simulation of a particle track 
 *  requires a particle mass, charge and momentum from the user. 
 *  The simulation of detector response smears the true intersections
 *  and creates TBHit objects for reconstruction. 
 * 
 *  Please note that the each gun object initializes an internal 
 *  random number generator to sample initial track state, 
 *  scatter angles and errors of position measurements. 
 *   
 *  Author: Benjamin Schwenker, Göttingen University 
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


/** Class TrackSnapshot
 * 
 * A TrackSnapshot is a pair (State/Surf) consiting of a 
 * ReferenceFrame (uvw tripod) and a vector (State) of 
 * local track parameters (du/dw,dv/dw,u,v) at the w=0 
 * surface. 
 *   
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class TrackSnapShot
{
 public:
  // Reference surface      
  ReferenceFrame Surf;
  // Track parameter vector
  CLHEP::HepMatrix State;    
};


class ParticleGun {

 public: 
   
  // Collimator X position in telescop coordinates (mm) 
  double positionX; 
  // Collimator Y position in telescop coordinates (mm)
  double positionY; 
  // Collimator Z position in telescop coordinates (mm)
  double positionZ; 
  // X Rotation of collimator plane (rad) 
  double rotationX; 
  // Y Rotation of collimator plane (rad) 
  double rotationY;
  // Z Rotation of collimator plane (rad)
  double rotationZ; 
   
  // Divergence of particle beam around beam axis in XZ plane (rad) 
  double divergenceX;
  // Divergence of particle beam around beam axis in YZ plane (rad) 
  double divergenceY;
  // Beam correlation coefficient for (dx/dz,x) pairs 
  double correlationX;
  // Beam correlation coefficient for (dy/dz,y) pairs 
  double correlationY;
  // Beam spot rms in x direction (mm)
  double spotsizeX;
  // Beam spot rms in y direction (mm)
  double spotsizeY;
  
  // Choose model of multiple scattering for the telescope ( Highland:0 or Moliere:1 )
  int mscmodel;   
  // Choose digitizer model for pixel sensors ( Mimosa26:0 or Gaussian:1 or Box:2)
  int digi_type;
  // Sensor hit detection efficiency [0,1]
  double efficiency;   

  
 public:
   
  // Constructor
  ParticleGun( );

  // Destructor 
  ~ParticleGun( );
         
  /** Simulate a charged particle trajectory
   *    
   *  Returns truth TBTrack object with truth track states at all crossed sensor along
   *  the trajectory.  
   */    
  TBTrack SimulateTrajectory(TBDetector& detector,double mass, double charge, double mom);
        
  /** Simulate hits along particle trajectory
   *    
   *  Smears truth intersections found in TruthTrack and fills 
   *  TBHit objects into the HitStore. HitStore should be empty. 
   */    
  void SimulateTrackHits(TBTrack& TruthTrack, std::vector<TBHit>& HitStore);
  
 private: 

  // Set seed for internal random number generator
  void SetSeed( int seed );

  // Internal random number generator (Mersenne Twister is safe choice)
  TRandom3 * myRng;
  
  /** Generate new particle 
   */     
  TrackSnapShot GenerateTrack(double charge, double mom);
    
  /** Simulate hit on detector
   *
   * The measured hit is just a random smearing of the true track intersection. 
   * The smearing can be done using different digitizer models(digi_type). 
   */  
  TBHit SimulateHit(CLHEP::HepMatrix& State, Det& adet, int digi_type);     
  
};


} // Namespace

#endif 

