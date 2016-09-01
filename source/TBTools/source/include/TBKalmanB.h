#ifndef TBKalmanB_H
#define TBKalmanB_H 1

// DEPFETTrackTools includes
#include "TBTrack.h"
#include "TBHit.h"
#include "GenericTrackModel.h"

// CLHEP includes 
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>


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
 * fit results and quality indicators. 
 * 
 * The track fit is startet by calling the function 'Fit' upon 
 * a TBTrack instance. The algorithm works as follows
 * 
 * 1) Compute all intersections of the reference trajectory in the 
 * TBTrack with all sub detectors (TE). Without a magnetic field, a
 * straight line track model is used for track extrapolation. 
 *                                                                              
 * 2) Perform a double digital filter (forward and backward pass) and 
 * compute optimal estimates at all surfaces. The current track
 * state is used to linearize the computation of transport matrices 
 * and material effects.           
 * 
 * 3) Smoothed track parameter estimates and covariances and chi2's are 
 * calculated for all sub detectors. Moreover, the reference trajectory 
 * is updated using the smoother output. 
 * 
 * 4) An outlier rejection can be performed using a cut on the smoothed 
 * hit chi2. In that case, one more iterations of the track fit (restart
 * at step 1) should be done.
 * 
 * As is, the tracking code makes a number of important assumptions:
 * 
 * a) No magnetic field in tracking; a straight line track model. 
 * b) Curvature q/p is not part of the state vector; user must supply a guess.   
 * c) Energy loss is not taken into account; 
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


//! Class FilterState
/* 
 * Class FilterState manages the filter states at a detector plane. 
 * 
 * The FilterState class stores [x,C] variable pairs for the 
 * estimated track states. The [x,C] pair is stored both for 
 * predicted and updated estimates. 
 * 
 * @Author B. Schwenker, University of Göttingen
 * <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */

class FilterState
{
 public:
  
  /* Predicted estimator 
   */
  CLHEP::HepMatrix Pr_x;     
  CLHEP::HepSymMatrix Pr_C; 
  
  /* Updated estimator 
   */
  CLHEP::HepMatrix Up_x;     
  CLHEP::HepSymMatrix Up_C;  
  
  /* Predicted Chi2, zero if no hit found 
   */  
  double Chi2Pred; 
};


/* Data type used to store the results of a filter pass. 
 */ 
typedef std::vector<FilterState> FilterStateVec;


class TBKalmanB {
  
 public:

  
  /** Default Constructor  
   */
  TBKalmanB(TBDetector& detector);

  /** Default d'tor  
   */
  ~TBKalmanB();

  /** Performs track fitting. Returns error flag. 
   */
  bool Fit(TBTrack& trk);

  /** Extrapolates the track seed to all planes. Returns error flag. 
   */
  bool ExtrapolateSeed(TBTrack& trk);
  
  /** Set number of iterations for filter
   */
  void SetNumIterations(int i){NumIt=i;}
  
  /** Set outlier removal cut 
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
  
  /** Run filter on track. Returns fit chi2. 
   *    
   *  Input: RefStateVec    : Vector of reference track parameters  
   *  Input: IDIR           : -1 for backward filter, 1 for forward filter 
   *  Output : RESULT       : Estimated [x,C] pairss for all detectors 
   */
  double FilterPass(TBTrack& trk, std::vector<int>& CrossedTEs, std::vector<CLHEP::HepMatrix>& RefStateVec ,int idir, FilterStateVec& Result); 
  
  /** Compute the weighted means of forward and backward filter estimated. In a 
   *  numerically robust way. Returns error flag.
   */
  bool GetSmoothedData( CLHEP::HepMatrix& xb, CLHEP::HepSymMatrix& Cb, CLHEP::HepMatrix& rf, CLHEP::HepSymMatrix& Cf,
                        CLHEP::HepMatrix& xs, CLHEP::HepSymMatrix& Cs);

  /** Set number of degrees of freedom in trk
   */
  void SetNdof(TBTrack& trk);
  
  int MAP_FORWARD(   double theta2, 
                        CLHEP::HepMatrix& xref, ReferenceFrame& Surf, 
                        CLHEP::HepMatrix& nxref, ReferenceFrame& nSurf, 
                        CLHEP::HepMatrix& x0, CLHEP::HepSymMatrix& C0
                     ); 
   
  int MAP_BACKWARD(   double theta2, 
                        CLHEP::HepMatrix& xref, ReferenceFrame& Surf, 
                        CLHEP::HepMatrix& nxref, ReferenceFrame& nSurf, 
                        CLHEP::HepMatrix& x0, CLHEP::HepSymMatrix& C0
                     ); 
  
   
 // Private Members -----------------
 private:
  
  // Iterate double filter
  int NumIt;
   
  // Project states to hit coord
  CLHEP::HepMatrix H; 
   
  // Outlier chisqu cut  
  double OutlierCut; 
  
  double mom; 
  double mass;
  double charge;
  
  GenericTrackModel* TrackModel; 

  int ndim; // dimension of state

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
