// ParticleGun implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Local includes
#include "ParticleGun.h"
#include "ThreeDModel.h"

#include "MaterialEffect.h"
#include "HelixTrackModel.h"
#include "StraightLineTrackModel.h"

// C++ includes
#include <iostream>
#include <cassert>
#include <cstdlib>

// ROOT includes
#include <TMath.h>

// CLHEP includes
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
using namespace CLHEP;

namespace depfet {

// Constructor
ParticleGun::ParticleGun( ) 
{
  // Provide default values for all 
  // gun settings 
  
  // X Rotation of collimator plane (rad) 
  rotationX = 0; 
  // Y Rotation of collimator plane (rad) 
  rotationY = 0; 
  // Z Rotation of collimator plane (rad)
  rotationZ = 0; 
  // Collimator X position in telescop coordinates (mm) 
  positionX = 0;
  // Collimator Y position in telescop coordinates (mm) 
  positionY = 0; 
  // Collimator Z position in telescop coordinates (mm) 
  positionZ = 0; 

  // RMS track slope in X (rad) 
  divergenceX = 0; 
  // RMS track slope in Y (rad) 
  divergenceY = 0; 
  // Correlation coefficient for (dx/dz,x) pairs 
  correlationX = 0;
  // Correlation coefficient for (dy/dz,y) pairs 
  correlationY = 0;
  // RMS of beam spot in x (mm)
  spotsizeX = 0;
   // RMS of beam spot in y (mm)
  spotsizeY = 0;  

  // Provide default value for trajectory 
  // propagations 
  
  //MSC model switch  (default:Highland)
  mscmodel=0;
  
  
  // Provide default value for sensor 
  // response simulation (->hits)
    
  // Choose digitizer model for pixel sensors ( Mimosa26:0 or Gaussian:1 or Box:2)
  digi_type=0;  
  
  // Sensor detection efficiency [0,1]
  efficiency = 1;
  // Random number generator - Mersenne Twister -> high quality random numbers 
  myRng = new TRandom3();
         
}
         
// Destructor
ParticleGun::~ParticleGun( ) {
  delete  myRng; 
} 

// Set random number generator seed
void ParticleGun::SetSeed(int seed) {
  myRng->SetSeed(seed);
}

/** Simulate a charged particle trajectory
 *    
 *  Returns truth TBTrack object with truth track states at all crossed sensor along
 *  the trajectory.  
 */    
TBTrack ParticleGun::SimulateTrajectory(TBDetector& detector, double mass, double charge, double mom)
{

  
  // A charged particle trajectory is a list of consecutive  
  // snapshots at all crossed detector planes.
  // 
  // Finally, the trajectory is mapped to a TBTrack object 
  // for convinience. 
  
  // Final particle trajectory object
     
  TBTrack TruthTrack(detector);
  TruthTrack.SetMass( mass );
  TruthTrack.SetCharge( charge );
  TruthTrack.SetMomentum( mom ); 
  
  // Init the track model for simulation
  GenericTrackModel* TrackModel; 
  
  HepVector field(3,0);
  field[0] = detector.GetBx();
  field[1] = detector.GetBy();
  field[2] = detector.GetBz(); 
  
  if ( field.norm() == 0 ) {
    TrackModel = new StraightLineTrackModel();    
  } else {
    TrackModel = new HelixTrackModel(field); 
    //TrackModel->m_kappa = charge/mom; 
  }     
  


  // Generate initial track snapshot at last collimator 
  // of the beam line.   
        
  TrackSnapShot MyTrack = GenerateTrack(charge, mom);  
  
  // Number of sub detectors
  int nsensor = detector.GetNSensors(); 
     
  // Count sensor planes in beam direction; Start at particle gun
  // sitting in front of the telescope at plane -1
  int ipl = -1;
   
  // Main tracking loop 
  do  { 
        
    //cout << " Track at plane " << ipl << " with pars " << MyTrack.State << endl;
    //cout << " Frame is " << endl; MyTrack.Surf.PrintHepMatrix();  
     
    // Scatter at a material or sensor plane 
     
    if (ipl >= 0 ) { 
      
      Det& current_det = detector.GetDet(ipl);
        
      // Store state at current detector  
      // -------------------------------
      TruthTrack.GetTE(ipl).GetState().SetPars(MyTrack.State); 
      TruthTrack.GetTE(ipl).SetCrossed(true);
       
      // Scatter on current detector
      // --------------------------- 
      double dudw_trk = MyTrack.State[0][0];
      double dvdw_trk = MyTrack.State[1][0];
      double u_trk = MyTrack.State[2][0]; 
      double v_trk = MyTrack.State[3][0];  
      
      // Traversed length in detector layer (mm)
      double l0 = current_det.GetThickness(u_trk,v_trk)*std::sqrt(1 + dudw_trk*dudw_trk + dvdw_trk*dvdw_trk);  
      // Radiation length of detector layer(mm) 
      double X0 = current_det.GetRadLength(u_trk,v_trk);
      
      // Scatter kinks (relative to flight direction of particle)      
      double kink_u = 0;
      double kink_v = 0;
      double theta2 = 0;
            
      if(  mscmodel==1  ) { 
        theta2 = materialeffect::GetScatterTheta2(l0, X0, mass, charge, mom );      
        kink_u = myRng->Gaus(0, TMath::Sqrt( theta2 ));
        kink_v = myRng->Gaus(0, TMath::Sqrt( theta2 ));       
      } else { // Highland model scattering
        theta2 = materialeffect::GetScatterTheta2(l0, X0, mass, charge, mom );      
        kink_u = myRng->Gaus(0, TMath::Sqrt( theta2 ));
        kink_v = myRng->Gaus(0, TMath::Sqrt( theta2 ));      
      }
       
      //cout << "det kinku " << kink_u << endl; 
      //cout << "det kinkv " << kink_v << endl; 
      
      // Scatter track ('in' state -> 'out' state)
      materialeffect::ScatterTrack(MyTrack.State, kink_u, kink_v);  
      
      // Also store scatter kinks in a TBHit object. 
      // Of course, this abuses the interface a bit.  
      TBHit MSCKink(current_det.GetDAQID(), kink_u, kink_v, theta2, theta2);
      TruthTrack.GetTE(ipl).SetHit(MSCKink); 
    }
         
    // Propagate track to next detector surface
    // ----------------------------------------
    
    int jpl = ipl+1;      
            
    // Exit condition   
    if (jpl == nsensor) {break;}
    
    // Next detector surface along beam line 
    ReferenceFrame next_surf = detector.GetDet(jpl).GetNominal();
               
    // Check that track really intersects surface
    if (! TrackModel->CheckHitsSurface(MyTrack.State, MyTrack.Surf, next_surf) ) {
      break;    
    } 
    
    // Now, compute the fligth length to next surface
    double length = TrackModel->GetSignedStepLength(MyTrack.State, MyTrack.Surf, next_surf);
    
    // Extrapolate half step 
    TrackModel->Extrapolate(MyTrack.State, MyTrack.Surf, 0.5*length); 
        
    // Scatter in air between detectors
    // ------------------------------- 
      
    double kink_u = 0;
    double kink_v = 0;
      
    if(mscmodel==1 ) { 
      // Sample the scattering angle in comoving frame      
      double theta2 = materialeffect::GetScatterTheta2(length, materialeffect::X0_air, mass, charge, mom ) ;   
      kink_u = myRng->Gaus(0, TMath::Sqrt( theta2 ));
      kink_v = myRng->Gaus(0, TMath::Sqrt( theta2 ));  
    } else { 
      // Sample the scattering angle in comoving frame      
      double theta2 = materialeffect::GetScatterTheta2(length, materialeffect::X0_air, mass, charge, mom ) ;   
      kink_u = myRng->Gaus(0, TMath::Sqrt( theta2 ));
      kink_v = myRng->Gaus(0, TMath::Sqrt( theta2 ));
    }
       
    //cout << "air kinku " << kink_u << endl; 
    //cout << "air kinkv " << kink_v << endl; 
        
    // Scatter track ('in' state -> 'out' state)
    materialeffect::ScatterTrack(MyTrack.State, kink_u, kink_v);    
          
    // Get state on next detector
    bool error = false; 
    HepMatrix next_state = TrackModel->Extrapolate(MyTrack.State, MyTrack.Surf, next_surf, error);  
    
    if (error) { // check for loopers
      break;      
    } 
    
    // Update track snapshot  
    ipl = jpl;  
    MyTrack.Surf = next_surf;
    MyTrack.State = next_state;  
    
               
  } while ( ipl < nsensor ); 
  
  delete TrackModel;  
 
  return TruthTrack;
}
        
/** Simulate hits along particle trajectory
 *    
 *  Smears truth intersections found in TruthTrack and fills 
 *  TBHit objects into the HitStore. HitStore should be empty. 
 */    
void ParticleGun::SimulateTrackHits(TBTrack& TruthTrack, std::vector<TBHit>& HitStore)
{
  
  // Number of detectors
  int nSensor = TruthTrack.GetNumTEs(); 
  
  // Force empty HitStore   
  HitStore.clear();
    
  for (int ipl = 0; ipl < nSensor; ++ipl) {
    
    // Skip sensors not crossed
    
    if ( !TruthTrack.GetTE(ipl).IsCrossed() ) continue; 
             
    // Get truth track intersect
    
    HepMatrix pars = TruthTrack.GetTE(ipl).GetState().GetPars();  
    double ut = pars[2][0]; 
    double vt = pars[3][0]; 
       
    // Check detector acceptance 
    
    Det& adet = TruthTrack.GetTE(ipl).GetDet();
    if ( ! adet.SensitiveCrossed(ut,vt) ) continue;
    
    // Check detector is sensitive 
    
    if ( adet.GetDeviceType() < 0 ) continue; 
    
    // Simulate detection efficiency
                  
    if ( myRng->Uniform() <= efficiency ) {
      TBHit hit = SimulateHit(pars, adet, digi_type);
      HitStore.push_back(hit); 
    }  
      
  }
  
  return;   
}
  
  
/** Generate new particle 
 */   
TrackSnapShot ParticleGun::GenerateTrack(double charge, double mom)
{
  
  // Particles are created in the plane of the particle
  // beam collimator opening 
   
  // First, we must construct of reference frame uvw plane for the
  // collimator opening 
  // The beam spot center is at u=v=0. The direction of the beam is 
  // the ew direction. 
  
  ReferenceFrame CollimatorFrame;
  
  // Position of center of collimator 
  HepVector CollimatorPosition(3);
  CollimatorPosition[0] = positionX; 
  CollimatorPosition[1] = positionY;
  CollimatorPosition[2] = positionZ;
  CollimatorFrame.SetPosition(CollimatorPosition); 
  
  // Rotation matrix of collimator plane  
  HepMatrix CollimatorRotation;
  double alpha = rotationX; 
  double beta = rotationY; 
  double gamma = rotationZ; 
  FillRotMatrixKarimaki(CollimatorRotation, alpha,beta,gamma);
  CollimatorFrame.SetRotation(CollimatorRotation);    
   
  // Next, we create a local track state in the collimator frame.
  HepMatrix GunState(5,1); 
  
  // Sample from a correlated particle position and direction 
  // variables 
  // 
  // Mean position and directions are always zero in collimator 
  // frame. 
  // 
  // Variables are decomposed into two uncorrelated pairs, namely 
  // (dx/dz,x) and (dy/dz,y). 
  //
  // We use a chelesky decomposition method to sample from the 
  // correlated pairs. 
   
  double X1 = myRng->Gaus(0, 1); 
  double X2 = myRng->Gaus(0, 1); 
  double Y1 = myRng->Gaus(0, 1); 
  double Y2 = myRng->Gaus(0, 1);

  double covX = correlationX*divergenceX*spotsizeX;
  double covY = correlationY*divergenceY*spotsizeY;  
  double DX = std::sqrt( std::pow(divergenceX,2)*std::pow(spotsizeX,2) - covX*covX );
  double DY = std::sqrt( std::pow(divergenceY,2)*std::pow(spotsizeY,2) - covY*covY );
  
  GunState[0][0] = divergenceX*X1;                     // dx/dz
  GunState[1][0] = divergenceY*Y1;                     // dy/dz
  GunState[2][0] = (covX*X1 + DX*X2) / divergenceX;    // x 
  GunState[3][0] = (covY*Y1 + DY*Y2) / divergenceY;    // y  
  GunState[4][0] = charge/mom;    // q/p
  
  // Finally, we put everything together :)
  TrackSnapShot MyTrack; 
  MyTrack.Surf = CollimatorFrame;
  MyTrack.State = GunState;
  
  return MyTrack;
}
    

    
/** Simulate hit on detector
 *
 * The measured hit is just a random smearing of the true track intersection. 
 */  
TBHit ParticleGun::SimulateHit(HepMatrix& State, Det& adet, int digi_type)
{
  // True particle intersection in local frame      
  double u = State[2][0];
  double v = State[3][0];  
  
  // Covariances of reconstructed hit
  double cov_u=0;
  double cov_v=0;   
   
  if ( digi_type == 2) {
    // Box smearing digitizer
    double pitch_u = adet.GetPitchU();
    double pitch_v = adet.GetPitchV();
    u += myRng->Uniform(-pitch_u/10., pitch_u/10.);
    v += myRng->Uniform(-pitch_v/10., pitch_v/10.);
    // Get column and row of seed pixel 
    int column=adet.GetColumnFromCoord(u,v);
    int row=adet.GetRowFromCoord(u,v);
    // Get the position of the center of the seed pixel
    u=adet.GetPixelCenterCoordU(row, column);
    v=adet.GetPixelCenterCoordV(row, column);
    cov_u = pow(pitch_u,2)/12.;  
    cov_v = pow(pitch_v,2)/12.;   
  } else if (digi_type == 1 ) { 
    // Gaussian smearing digitizer 
    double reso_u = adet.GetResolutionU(); 
    double reso_v = adet.GetResolutionV(); 
    u += myRng->Gaus(0, reso_u); 
    v += myRng->Gaus(0, reso_v); 
    cov_u = pow(reso_u,2); 
    cov_v = pow(reso_v,2); 
  } else {
    // Binary Mimosa26 digitizer parameters
    double pitch_u = adet.GetPitchU();
    double pitch_v = adet.GetPitchV();
    double rsmear = 0.0; 
    double borderU = 0.5*pitch_u;
    double borderV = 0.5*pitch_v;     
    
    // Start binary Mimosa26 digitizer
    u += myRng->Gaus( 0, rsmear );
    v += myRng->Gaus( 0, rsmear ); 
    // Get column and row of seed pixel 
    int column=adet.GetColumnFromCoord(u,v);
    int row=adet.GetRowFromCoord(u,v);
    // Get the position of the center of the seed pixel
    double u_center=adet.GetPixelCenterCoordU(row, column);
    double v_center=adet.GetPixelCenterCoordV(row, column);
    double deltaU = u - u_center; 
    double deltaV = v - v_center;  
     
    if ( deltaU < -borderU/2. ) {
      u = u_center - pitch_u/2.;
      cov_u = pitch_u*pitch_u/48 + rsmear*rsmear; 
    } else if ( deltaU > borderU/2. ) {
      u = u_center + pitch_u/2.;
      cov_u = pitch_u*pitch_u/48 + rsmear*rsmear; 
    } else {
      u = u_center;
      cov_u = pitch_u*pitch_u/48 + rsmear*rsmear; 
    } 
    
    if ( deltaV < -borderV/2. ) {
      v = v_center - pitch_v/2.;
      cov_v = pitch_v*pitch_v/48 + rsmear*rsmear; 
    } else if ( deltaV > borderV/2. ) {
      v = v_center + pitch_v/2.;
      cov_v = pitch_v*pitch_v/48 + rsmear*rsmear; 
    } else {
      v = v_center;
      cov_v = pitch_v*pitch_v/48 + rsmear*rsmear; 
    } 
        
  }
   
  return TBHit(adet.GetDAQID(), u, v, cov_u, cov_v) ; 
}   

          
} // Namespace;
