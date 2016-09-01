#ifndef TBKalmanSRIF_H
#define TBKalmanSRIF_H 1

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TBHit.h"
#include "GenericTrackModel.h"

// CLHEP includes 
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

namespace depfet {

 
//! Class TBKalmanSRIF 
/* 
 * The TBKalmanSRIF class performs fitting of charged particle tracks
 * in a pixel tracking telescope. The track fit is based on a  
 * Kalman Smoother and provides optimal track parameters at all sub 
 * detectors. The track fit takes into account multiple scattering 
 * in sub detectors and air gaps between the sensors. Moreover, the 
 * track fit can deal with arbitrary 3D orientations of sub detectors.  
 * 
 * The TBKalmanSRIF class operates on TBTrack objects. Each instance 
 * of TBTrack holds all required input data to carry out the fit. 
 * After the fit is completed, the TBTrack instance also holds 
 * fit results and quality indicators. 
 * 
 * The track fit is startet by calling the function 'Fit' upon 
 * a TBTrack instance. The algorithm works as follows
 * 
 * 1) Compute all intersections of the reference trajectory in the 
 * TBTrack with all sub detectors. Without a magnetic field, a
 * straight line track model is used for track extrapolation. 
 *                                                                              
 * 2) Perform a double Kalman filter (forward and backward pass) and 
 * compute optimal estimates at all surfaces. The track fit is linearized 
 * the around the reference trajectory. It means that transport matrices 
 * and material effects are calculated for the reference trajectory.            
 * 
 * 3) Smoothed track parameter estimates and covariances and chi2's are 
 * calculated for all sub detectors. Moreover, the reference trajectory 
 * is updated using the smoother output. 
 * 
 * 4) An outlier rejection can be performed using a cut on the smoothed 
 * hit chi2. In that case, one more iterations of the track fit (restart
 * at step 1) should be done.
 * 
 * As is, the tracking code makes a number of assumptions:
 * 
 * a) The tracking detector is a linear array of pixel sensors
 * b) No magnetic field in tracking; a straight line track model. 
 * c) The user must supply a momentum measurment (beam energy).   
 * d) No Bethe-Bloch energy loss and no Bremsstrahlung energy loss.  
 * 
 * The last assumption is save if the relative total energy loss dE/E 
 * is small (<1%). Typically, this requires beam energies (>>100MeV) 
 * and a small material budget (~1% X0) per sub detector. 
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
 * The KalFilterDet class stores [R,z] variable pairs for the 
 * estimated track states. If [C,x] are estimated covariance 
 * and track parameters, the relations are   
 *  
 * a) Y=R^t*R , R is a square root of the information matrix Y=1/C
 * b) z=Rx  
 * 
 * The pair [R,z] is the natural building block for a Square Root 
 * Information Filter (SRIF) mechanization of Kalman's algorithm.   
 *  
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class KalFilterDet
{
 public:
  
  /* Predicted estimator 
   */
  CLHEP::HepMatrix Pr_z;     // 4x1 matrix
  CLHEP::HepMatrix Pr_R;     // 4x4 matrix
  
  /* Updated estimator 
   */
  CLHEP::HepMatrix Up_z;     // 4x1 matrix
  CLHEP::HepMatrix Up_R;     // 4x4 matrix
  
  /* Predicted Chi2, zero if no hit found 
   */  
  double Chi2Pred; 
};


/* Data type used to store the results of a filter pass. 
 */ 
typedef std::vector<KalFilterDet> KalFilterVec;

/* Data type used to store local track parameters along reference 
 * trajectory. 
 */ 
typedef std::vector<CLHEP::HepMatrix> REFTrack;

class TBKalmanSRIF {
  
 public:
  
  /** Default Constructor  
   */
  TBKalmanSRIF();

  /** Default Destructor 
   */
  ~TBKalmanSRIF();
  
  /** Performs track fitting. Returns error flag. 
   */
  bool Fit(TBTrack& trk);
  
  /** Extrapolates the track seed to all planes. Returns error flag. 
   */
  bool ExtrapolateSeed(TBTrack& trk);
  
  /** Set number of iterations for Kalman filter
   */
  void SetNumIterations(int i){NumIt=i;}
  
  /** Set outlier removal cut (Smoother ChiSqu)
   */
  void SetOutlierCut(double cut){OutlierCut=cut;}
  
  /** Set outlier removal cut (Smoother ChiSqu)
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
  double GetChi2Increment(CLHEP::HepMatrix& p, CLHEP::HepSymMatrix& C, TBHit& hit); 
  
  /** Returns the chi2 increment
   */
  double GetChi2Increment( CLHEP::HepMatrix& r, CLHEP::HepMatrix& H, CLHEP::HepSymMatrix& C, CLHEP::HepSymMatrix& V); 
  
  /** Returns the predicted chi2 for hit
   */
  double GetPredictedChi2(CLHEP::HepMatrix& p, CLHEP::HepSymMatrix& C, TBHit& hit); 
  
  /** Returns the predicted chi2
   */
  double GetPredictedChi2( CLHEP::HepMatrix& r, CLHEP::HepMatrix& H, CLHEP::HepSymMatrix& C, CLHEP::HepSymMatrix& V); 
    
  CLHEP::HepMatrix& GetHMatrix()
    { return H;} 
    
 // Private Methods -----------------
 private:
    
  /** Run Kalman filter on track. Returns fit chi2. 
   *    
   *  Input: StateVec    : Vector of reference track parameters  
   *  Input: IDIR        : -1 for backward filter, 1 for forward filter 
   *  Output : RESULT    : Estimated [R,z] pairs for all detectors 
   */
  double KalmanFilterPass(TBTrack& trk, std::vector<int>& CrossedTEs, REFTrack& StateVec, int idir, KalFilterVec& Result); 
   
  /** Compute the weighted means of forward and backward filter estimated. In a 
   *  numerically robust square root way. Returns error flag.
   */
  bool GetSmoothedData( CLHEP::HepMatrix& zb, CLHEP::HepMatrix& Rb, CLHEP::HepMatrix& zf, CLHEP::HepMatrix& Rf,
                        CLHEP::HepMatrix& Smoothed_State, CLHEP::HepSymMatrix& Smoothed_Cov);

  /** Compute a priori predicted estimate on first sensor. Returns reference state at first sensor .
   */
  CLHEP::HepMatrix ComputeBeamConstraint( CLHEP::HepMatrix& z0, CLHEP::HepMatrix& R0, ReferenceFrame& FirstSensorFrame, double mass, double mom, double charge);
  
  /** Returns a proper start column for Householder triangularization 
   *  
   * The method tries to triangularize first maxcol columns of A 
   * by row permutations. It returns the first column, where this 
   * is not enough.   
  */
  int StartupHouseholder(CLHEP::HepMatrix& A, int maxcol); 
  
  int MAP_FORWARD(   double theta2, 
                        CLHEP::HepMatrix& xref, ReferenceFrame& Surf, 
                        CLHEP::HepMatrix& nxref, ReferenceFrame& nSurf, 
                        CLHEP::HepMatrix& z0, CLHEP::HepMatrix& R0
                     ); 
   
  int MAP_BACKWARD(   double theta2, 
                        CLHEP::HepMatrix& xref, ReferenceFrame& Surf, 
                        CLHEP::HepMatrix& nxref, ReferenceFrame& nSurf, 
                        CLHEP::HepMatrix& z0, CLHEP::HepMatrix& R0
                     ); 
  
 // Private Members -----------------
 private:
  
  // Iterate double Kalman filter
  int NumIt;
  
  // Outlier chisqu cut  
  double OutlierCut; 
  
  // Project states to hit coord
  CLHEP::HepMatrix H;  
  
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
