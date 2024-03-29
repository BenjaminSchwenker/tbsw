#ifndef TBKalmanMSC_H
#define TBKalmanMSC_H 1

// TBTools includes
#include "TBTrack.h"
#include "TBHit.h"
#include "GenericTrackModel.h"

//other includes includes
#include <TMath.h>


namespace depfet {

 
//! Class TBKalmanMSC 
/* 
 * The TBKalmanMSC class performs fitting of charged particle tracks
 * in a pixel tracking telescope. The track fit is based on a Kalman
 * Filter and takes into account multiple scattering in sensors and 
 * air gaps between sensors. The fitter performs either a forward or
 * a backward pass of the Kalman filter. No smoothing of track states 
 * on inner sensors is done.   
 *
 * The TBKalmanMSC class operates on TBTrack objects. Each instance 
 * of TBTrack holds all required input data to carry out the fit. 
 * After the fit is completed, the TBTrack instance also holds 
 * fit results and quality indicators. 
 * 
 * The track fit is startet by calling the function 'ProcessTrack'
 * upon a TBTrack instance. The algorithm works as follows:
 * 
 * 1) Compute all intersections of the reference trajectory in the 
 * TBTrack with all sub detectors (TE). Without a magnetic field, a
 * straight line track model is used for track extrapolation. 
 *                                                                              
 * 2) Perform a forward or backward  Kalman filter pass and compute 
 * predicted and updated state estimates at all detectors. 
 * The track fit is linearized the around the reference trajectory. 
 * It means that transport matrices and material effects are calculated 
 * for the reference trajectory.            
 * 
 * 3) Either biased or unbiased states are copied to the track elements 
 * of the TBTrack object.  
 * 
 * As is, the tracking code makes a number of important assumptions:
 * 
 * a) No magnetic field in tracking; a straight line track model. 
 * b) Curvature q/p is not part of the state vector; user must supply guess.   
 * c) Energy loss is not taken into account; minor effect w/o magnetic.
 * d) Only pixel measurments used. So far, i never tried strips. 
 * 
 * These assumptions are save as long as there is no magnetic field and 
 * the relative total energy loss dE/E is small (<1%). Typically, this  
 * requires beam energies (>100MeV) and a small material budget (<1%X0)
 * per sub detector. Strip detectors could be integrated using appropriatly 
 * large errors for the v axis. 
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

//! Class KalFilterDet
/* 
 * Class KalFilterDet manages the results of a Kalman filter pass 
 * for a particular sub detector. The class can be symmetrically
 * used for the forward and backward filter pass. The detectors
 * are uniquily identified by their plane number, i.e. position 
 * along the beam line. 
 * 
 * The KalFilterDet class stores [x,C] variable pairs for the 
 * estimated track states. The [x,C] pair is stored both for 
 * predicted and updated KF estimates. 
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class KalFilterDetMSC
{
 public:
  KalFilterDetMSC():Pr_x(TrackState::Zero()),Pr_C(TrackStateCovariance::Zero()),
    Up_x(TrackState::Zero()),Up_C(TrackStateCovariance::Zero()){ }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /* Predicted estimator 
   */
  TrackState Pr_x;    
  TrackStateCovariance Pr_C;     
  
  /* Updated estimator 
   */
  TrackState Up_x;     
  TrackStateCovariance Up_C;    
  
  /* Predicted Chi2, zero if no hit found 
   */  
  double Chi2Pred; 
};

/* Data type used to store the results of a filter pass. 
 */ 
typedef std::vector<KalFilterDetMSC> KalFilterVecMSC;

/* Data type used to store local track parameters along reference 
 * trajectory. 
 */ 
typedef std::vector<TrackState> REFTrack;

typedef Eigen::Matrix<double,2,5> StateHitProjector;

class TBKalmanMSC {
  
 public:
  
  /** Default Constructor  
   */
  TBKalmanMSC(TBDetector& detector);

  /** Default Destructor 
   */
  ~TBKalmanMSC();
  
  /** Performs track fitting. Returns track chisqu. 
   */
  bool ProcessTrack(TBTrack& trk, int dir, bool biased);
  
  /** Get scatter kink angles (as 2x1 matrix) 
   */
  TrackScatterKinks GetScatterKinks(TBTrackState& InState, TBTrackState& OutState);

  /** Get covariance of scatter kinks (as 2x2 sym matrix)
   */
  TrackScatterKinksCovariance GetScatterKinkCov(TBTrackState& InState, TBTrackState& OutState);

  /** Help function to get thetas from slope vector
   */
  TrackScatterKinks slopestotheta(double slopes[4]);
  
  const StateHitProjector& GetHMatrix()
    { return H;}


 // Private Methods -----------------
 private:
  
  /** Run Kalman filter on track. Returns fit chi2. 
   *  
   *  Input: IDIR : -1 for backward filter, 1 for forward filter
   *  ISTART   : Plane number of start detector 
   *  ISTOP    : Plane number of last detector 
   *  Output : RESULT : Estimated [R,z] pairs for all detectors
   *  RSTATEVEC : Local states along reference trajectory 
   */
  double KalmanFilter(TBTrack& trk, int idir, int istart, int istop, KalFilterVecMSC& Result, REFTrack& RStateVec); 

  /** Set number of degrees of freedom in trk
   */
  void SetNdof(TBTrack& trk);
   
 // Private Members -----------------
 private:
  
  // Project states to hit coord
  StateHitProjector H;

  // Initial track state vector to start fitting 
  TrackState x0;

  // Initial track state covariacne matric to start fitting
  TrackStateCovariance C0;

  GenericTrackModel* TrackModel;

  int ndim; // dimension of state
 
  double mass;
  double charge;
  

};
 

} // Namespace

#endif 
