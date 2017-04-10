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
 * TBTrack with all sub detectors (TE).
 *                                                                              
 * 2) Perform  double filter pass (forward and backward) and 
 * compute optimal predicted estimates at all sensor planes. The current track
 * state is used to linearize the computation of transport matrices 
 * and material effects.           
 * 
 * 3) Smoothed track parameter estimates and covariances and chi2's are 
 * calculated for all sensor planes. Moreover, the reference trajectory 
 * is updated using the smoother output. 
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
 * predicted estimate only 
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
  
  /** Propagte state from trackelement te to next track element nte 
   */
  int PropagateState(TBTrackElement& te, TBTrackElement& nte, CLHEP::HepMatrix& xref, CLHEP::HepMatrix& nxref, CLHEP::HepMatrix& x0, CLHEP::HepSymMatrix& C0);
  
  /** Filters a new hit. Returns predicted chi2 
   */
  double FilterHit(TBHit& hit, CLHEP::HepMatrix& xref, CLHEP::HepMatrix& x0, CLHEP::HepSymMatrix& C0);  

  /** Get underlying track model
   */
  GenericTrackModel* GetTrackModel() {return TrackModel;}  
   
  CLHEP::HepMatrix& GetHMatrix()
    { return H;} 

  /** Set number of degrees of freedom in trk
   */
  void SetNdof(TBTrack& trk);
     
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
