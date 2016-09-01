// Implementation of TrackInputProvider class.
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// DEPFETTrackTools includes
#include "TrackInputProvider.h"

// LCIO includes
#include <lcio.h>
#include <IMPL/TrackerHitImpl.h>
	
using namespace lcio;	
using namespace CLHEP;
using namespace std;

namespace depfet {
 	
/** Build LCIO Track 
 * 
 * The lcio::Track class is used to store reco tracks persistently in 
 * LCIO files.  
 */
TrackImpl * TrackInputProvider::MakeLCIOTrack( TBTrack& trk )
{ 
  // Create a new lcio track object 
  TrackImpl * fittrack = new TrackImpl();
  
  // Copy particle hypothesis   
  float refpoint[3] = { trk.GetMass(), trk.GetCharge(), trk.GetMomentum() };
  fittrack->setReferencePoint(refpoint);
  
  // Copy quality indicators
  fittrack->setChi2(trk.GetChiSqu());   
  
  // Reference plane number 
  fittrack->setNdf( trk.GetReferenceState().GetPlane()  );           
  
  // Copy reference track 
  HepMatrix refstate = trk.GetReferenceState().GetPars();
   
  fittrack->setTanLambda(refstate[0][0]);    // This is u'
  fittrack->setPhi(refstate[1][0]);          // This is v' 
  fittrack->setD0(refstate[2][0]);           // This is u 
  fittrack->setZ0(refstate[3][0]);           // This is v                                   
  fittrack->setOmega(refstate[4][0]);        // This is q/p
   
  // Store hits in lcio track  
  int nTE = trk.GetNumTEs(); 
  
  for(int iTE=0; iTE<nTE; ++iTE) {
    if ( trk.GetTE(iTE).HasHit() ) {
      TrackerHit *lciohit = trk.GetTE(iTE).GetHit().GetRawHit();    
      fittrack->addHit(lciohit);  
    }
  }
  
  return fittrack; 
} 

/** Build TBTrack from lcio::Track
 * 
 * The lcio::Track class is used to store reco tracks persistently in 
 * LCIO files.  
 * 
 */
TBTrack TrackInputProvider::MakeTBTrack( lcio::Track * lciotrk, TBDetector& detector )
{  
  
  // Create track  
  TBTrack trk(detector);
  
  // Set particle hypothesis   
  const float* refpoint = lciotrk->getReferencePoint();
  trk.SetMass(refpoint[0]);
  trk.SetCharge(refpoint[1]);
  trk.SetMomentum(refpoint[2]);  
       
  // Set track chisqu  
  trk.SetChiSqu(lciotrk->getChi2());
  
  // Create reference track state (seed)
  HepMatrix Pars(5,1,0); 
  Pars[0][0] = lciotrk->getTanLambda();    // This is u'
  Pars[1][0] = lciotrk->getPhi();          // This is v' 
  Pars[2][0] = lciotrk->getD0();           // This is u 
  Pars[3][0] = lciotrk->getZ0();           // This is v
  Pars[4][0] = lciotrk->getOmega();        // This is q/p
  
  TBTrackState Seed;
  Seed.Pars = Pars;
  Seed.SetPlane( lciotrk->getNdf() ); 
  trk.SetReferenceState(Seed);
    
  // Add hits to track  
  std::vector<TrackerHit*>  trackhits = lciotrk->getTrackerHits();      
  int nHit = (int) trackhits.size();
  
  for(int ihit=0;ihit<nHit;ihit++)  { 
    // Create a TBHit     
    TrackerHitImpl * lciohit = dynamic_cast<TrackerHitImpl*>(trackhits.at(ihit));  
    TBHit Hit(lciohit);   
    // Add TBHit to TBTrack  
    int daqid = Hit.GetDAQID();
    int ipl = detector.GetPlaneNumber(daqid);
    trk.GetTE(ipl).SetHit(Hit); 
  } 
   
  return trk; 
}  
 	

} // Namespace;
