#ifndef KalmanAlignmentAlgorithm2_H
#define KalmanAlignmentAlgorithm2_H 1

// TBTools includes
#include "TBDetector.h"
#include "AlignableDet.h"

#include "GenericTrackModel.h"

// ROOT includes
#include "TFile.h"

namespace depfet {



//! Class KalmanAlignmentAlgorithm2 
/*  
 *  The KalmanAlignmentAlgorithm2 class fits the alignment constants of 
 *  HEP detectors using track residuals. The method is iterative, 
 *  based on the Kalman filter equations. The Kalman filter equations 
 *  offer the possibility to use prior information about the alignment 
 *  from mechanical measurements or previous alignment steps. And they 
 *  easily allow to fix alignment parameters to define a global frame 
 *  of reference. The method is also suitable for alignment of sub 
 *  detector relative to a different detector. 
 *    
 *  The implementation of the Kalman Alignment Algorithm (KAA) follows 
 *  Rudi Fruehwirth's paper "Estimation of Alignment Parameters, using 
 *  Kalman Filter with annealing", CMS Note 2002/008. 
 *   
 *  @Author B. Schwenker, University of Göttingen
 *  <mailto:benjamin.schwenker@phys.uni-goettingen.de>
 */
  
class KalmanAlignmentAlgorithm2 {
  
 public:
  
  /** Default Constructor  
   */
  KalmanAlignmentAlgorithm2();
    
  /** Returns new reference frame after applying sensor alignment parameters. 
   */
  static ReferenceFrame ApplyAlignment(const ReferenceFrame& nominal, const SensorAlignmentParameters& alignPars); 
   
  /** Align detector with align constants. 
   */
  static void AlignDetector(TBDetector& detector, const  AlignableDet & alignconst);
  
  /** Performs alignment fit. Returns alignment results. 
   */
  AlignableDet Fit(TBDetector& detector, TFile * AlignmentData, AlignableDet initialAlignState, 
                   int maxTracks, int annealingTracks, double annealingFactor, double  pValueCut, double deviationCut, bool useBC, int logLevel);
  
 private:
  
  /** Jacobian matrix, derivatives of measurement equation f=(fu,fv) to 
   *  the six alignment parameters ar=(dx, dy, dz, dalpha, dbeta, dgamma).
   *
   *  The Jacobian matrix is evaluated for track parematers p0=(tu, tv, u, v)
   *  at given sensor position and rotation.   
   */
  SensorAlignmentJacobian Jacobian_Alignment(const TrackState & p0, const Matrix3d & Rot ) const;
   
};
 

} // Namespace

#endif 
