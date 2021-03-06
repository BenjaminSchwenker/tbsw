#ifndef TBKalmanB_H
#define TBKalmanB_H 1

// TBTools includes
#include "TBTrack.h"
#include "TBHit.h"
#include "GenericTrackModel.h"

namespace depfet {
//! Class KalFilterDet
/*
 * Class KalFilterDet contains a Kalman filter estimate of 
 * a local track state. It contains the state parameter 
 * vector and its covariance matrix.
 *
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class KalFilterDet
{
 public:
  KalFilterDet():Pr_x(TrackState::Zero()),Pr_C(TrackStateCovariance::Zero()){ }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /* Predicted estimator
   */
  TrackState Pr_x;
  TrackStateCovariance Pr_C;

};
}
#include<Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(depfet::KalFilterDet)

namespace depfet {

 
//! Class TBKalmanB 
/* 
 * The TBKalmanB class performs fitting of charged particle tracks
 * in a pixel tracking telescope. The track fit is based on recursive
 * Bayesian filtering. It takes into account multiple 
 * scattering both in detector modules and air gaps between detectors. 
 * 
 * The TBKalmanB class operates on TBTrack objects. Each instance 
 * of TBTrack holds all required input data to carry out the fit.  
 * After the fit is completed, the TBTrack instance also holds 
 * fit results and fit quality indicators. 
 * 
 * The track fit is startet by calling the function 'Fit' upon 
 * a TBTrack instance. The algorithm works as follows
 * 
 * 1) Compute all intersections of the reference trajectory in the 
 * TBTrack with all sub detectors (TE).
 *                                                                              
 * 2) Perform  double filter pass (forward and backward) and 
 * compute optimal predicted estimates at all sensor planes. The 
 * local reference state is used to linearize the computation of 
 * transport matrices and material effects.           
 * 
 * 3) Smoothed track parameter estimates and covariances and chi2's are 
 * calculated for all sensor planes. Moreover, the reference trajectory 
 * is updated using the smoother output. 
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */


/* Data type used to store the results of a filter pass. 
 */ 
typedef std::vector<KalFilterDet> FilterStateVec;
typedef Eigen::Matrix<double,2,5> StateHitProjector;

class TBKalmanB {
  
 public:

  
  /** Default Constructor  
   */
  TBKalmanB(TBDetector& detector);
  
  /** Default Destructor 
   */
  ~TBKalmanB();
  
  /** Performs track fitting. Returns error flag. 
   * 
   * The smoothPlane variable control the smoother pass. Use a value of -1 for smoothing all track elements.
   * Alternativly, the user can put the number (>=0) of one specific track element to restrict the smoothing
   * to this track element only. In this case, the track state and covariance matrix of the other track elements will 
   * be undefined.    
   */
  bool Fit(TBTrack& trk, int smoothPlane);

  /** Extrapolates the track seed to all planes. Returns error flag. 
   */
  bool ExtrapolateSeed(TBTrack& trk);
  
  /** Set number of iterations for filter
   */
  void SetNumIterations(int i){NumIt=i;}
  
  /** Set use beam constraint
   */
  void SetUseBeamConstraint(bool useBC){m_useBC=useBC;}
  
  /** Set constraint on beam divergence in rad
   */
  void SetBeamDivergenceX(double div){m_divx=div;}

  /** Set constraint on beam divergence in rad
   */
  void SetBeamDivergenceY(double div){m_divy=div;}

  /** Set constraint on beam size in mm
   */
  void SetBeamSizeX(double size){m_sizex=size;}

  /** Set constraint on beam size in mm
   */
  void SetBeamSizeY(double size){m_sizey=size;}

  /** Set constraint on beam coorelation betwween x and dx/dz 
   */
  void SetBeamCorrelationX(double corr){m_corrx=corr;}

  /** Set constraint on beam coorelation betwween y and dy/dz 
   */
  void SetBeamCorrelationY(double corr){m_corry=corr;}
  
  /** Returns the chi2 increment
   */
  double GetChi2Increment(const TrackState& p, const TrackStateCovariance& C, const TBHit& hit);
  
  /** Returns the chi2 increment
   */
  double GetChi2Increment(const Vector2d& r, const StateHitProjector& H, const TrackStateCovariance& C, const Matrix2d& V);
    
  /** Returns the predicted chi2 for hit
   */
  double GetPredictedChi2(const TrackState& p, const TrackStateCovariance& C, const TBHit& hit);
  
  /** Returns the predicted chi2
   */
  double GetPredictedChi2(const Vector2d& r, const StateHitProjector& H, const TrackStateCovariance& C, const Matrix2d& V);
  
  /** Propagte state from trackelement te to next track element nte 
   */
  int PropagateState(TBTrackElement& te, TBTrackElement& nte, TrackState& xref, TrackState& nxref, TrackState& x0, TrackStateCovariance& C0);
  
  /** Filters a new hit. Returns predicted chi2 
   */
  double FilterHit(const TBHit& hit, const TrackState& xref, TrackState& x0, TrackStateCovariance& C0);  

  /** Get underlying track model
   */
  GenericTrackModel* GetTrackModel() {return TrackModel;}  
   
  const StateHitProjector& GetHMatrix()
    { return H;} 

  /** Set number of degrees of freedom in trk
   */
  void SetNdof(TBTrack& trk);
     
 // Private Methods -----------------
 private:
  
  /** Compute a priori predicted estimate on first sensor. Returns reference state at first sensor.
   */
  TrackState ComputeBeamConstraint( ReferenceFrame& FirstSensorFrame, double mom, double charge);
  
  /** Run filter on track. Returns fit chi2. Failure if chi2<0.
   *    
   *  Input: trk            : Track with all associated hits and seed parameters
   *  Input: CrossedTEs     : Vector with the plane numbers of intersected senors planes
   *  Input: RefStateVec    : Vector of reference track parameters  
   *  Input: IDIR           : -1 for backward filter, 1 for forward filter 
   *  Output : RESULT       : Estimated [x,C] pairs for all detectors 
   */
  double FilterPass(TBTrack& trk, std::vector<int>& CrossedTEs, std::vector<TrackState>& RefStateVec ,int idir, FilterStateVec& Result); 
  
  /** Compute the weighted means of forward and backward filter estimated. In a 
   *  numerically robust way. Returns error flag.
   */
  bool GetSmoothedData( const TrackState& xb, const TrackStateCovariance& Cb, const TrackState& rf, const TrackStateCovariance& Cf,
                        TrackState& xs, TrackStateCovariance& Cs);
  
  int MAP_FORWARD(double theta2,
                        const TrackState& xref, const ReferenceFrame& Surf, const ReferenceFrame& nSurf,
                        TrackState& x0, TrackStateCovariance& C0
                     );
   
  int MAP_BACKWARD(double theta2,
                        const TrackState &xref, const ReferenceFrame &Surf,
                        const TrackState &nxref, const ReferenceFrame &nSurf,
                        TrackState& x0, TrackStateCovariance& C0
                     );
  
   
 // Private Members -----------------
 private:
  
  // Iterate double filter
  int NumIt;
   
  // Project states to hit coord
  StateHitProjector H;
   
  // Initial track state vector to start fitting 
  TrackState x0_init;

  // Initial track state covariacne matric to start fitting
  TrackStateCovariance C0_init;

  double mass;
  double charge;
  
  GenericTrackModel* TrackModel; 
  
  // Use beam constraint in fitter
  bool m_useBC; 
  double m_divx;      // divergence, rad
  double m_divy;      // divergence, rad
  double m_sizex;     // spot size, mm
  double m_sizey;     // spot size, mm 
  double m_corrx;     // beam correlation coeff  
  double m_corry;     // beam correlation coeff 
        
};
 

} // Namespace

#endif 
