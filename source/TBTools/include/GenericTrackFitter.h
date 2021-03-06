#ifndef GenericTrackFitter_H
#define GenericTrackFitter_H 1

// TBTools includes
#include "TBTrack.h"
#include "TBHit.h"
#include "TBKalmanB.h"

namespace depfet {

 
//! Class GenericTrackFitter
/* 
 * The GenericTrackFitter class defines an user interface for track 
 * fitting in the TBSW software. The main purpose of this class is 
 * to provide an easy way to change the track fitting method used 
 * for track reconstruction. 
 *
 * The GenericTrackFitter class operates on TBTrack objects. The 
 * class TBTrack holds all required input data to perform the fit. 
 * After the fit is completed, the TBTrack instance also holds the
 * fit results and quality indicators.  
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
  

class GenericTrackFitter {
  
 public:

  /** Default Constructor  
   */
  GenericTrackFitter(TBDetector& detector) { m_TrackFitter = new TBKalmanB(detector); }
 
  ~GenericTrackFitter() { delete m_TrackFitter; }

  /** Performs track fitting. Returns error flag. 
   *   
   *  The smoothPlane variable control the smoother pass. Use a value of -1 for smoothing all track elements.
   *  Alternativly, the user can put the number (>=0) of one specific track element to restrict the smoothing
   *  to this track element only. In this case, the track state and covariance matrix of the other track elements will 
   *  be undefined.    
   */
  bool Fit(TBTrack& trk, int smoothPlane=-1){return m_TrackFitter->Fit(trk, smoothPlane);}
  
  /** Extrapolates the track seed to all planes. Returns error flag. 
   */
  bool ExtrapolateSeed(TBTrack& trk){return m_TrackFitter->ExtrapolateSeed(trk);} 
  
  /** Set number of iterations for Kalman filter
   */
  void SetNumIterations(int i){m_TrackFitter->SetNumIterations(i);}
   
  /** Set use beam constraint
   */
  void SetUseBeamConstraint(bool useBC){m_TrackFitter->SetUseBeamConstraint(useBC);}
  
  /** Set constraint on beam divergence in rad
   */
  void SetBeamDivergenceX(double div){m_TrackFitter->SetBeamDivergenceX(div);}

  /** Set constraint on beam divergence in rad
   */
  void SetBeamDivergenceY(double div){m_TrackFitter->SetBeamDivergenceY(div);}

  /** Set constraint on beam size in mm
   */
  void SetBeamSizeX(double size){m_TrackFitter->SetBeamSizeX(size);}

  /** Set constraint on beam size in mm
   */
  void SetBeamSizeY(double size){m_TrackFitter->SetBeamSizeY(size);}

  /** Set constraint on beam coorelation betwween x and dx/dz 
   */
  void SetBeamCorrelationX(double corr){m_TrackFitter->SetBeamCorrelationX(corr);}

  /** Set constraint on beam coorelation betwween y and dy/dz 
   */
  void SetBeamCorrelationY(double corr){m_TrackFitter->SetBeamCorrelationY(corr);}

  /** Get the projector matrix to exract the position from track state
   */ 
  const StateHitProjector& GetHMatrix(){return  m_TrackFitter->GetHMatrix();} 


  /** Returns the predicted local chi2 for the given hit and track parameter estimate
   */
  double GetPredictedChi2(const TrackState& p, const TrackStateCovariance& C, const TBHit& hit) {return m_TrackFitter->GetPredictedChi2(p,C,hit); } 
  
 private:
  
  // This is the switch to choose the fitter
  //TBKalmanSRIF * m_TrackFitter;  
  TBKalmanB * m_TrackFitter;  
  
};
 
} // Namespace

#endif 
