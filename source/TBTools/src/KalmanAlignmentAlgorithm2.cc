// KalmanAlignmentAlgorithm2 implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// C++ includes 
#include <iostream>

// TBTool includes
#include "KalmanAlignmentAlgorithm2.h"
#include "AlignEvent.h"
#include "KalmanAlignmentInputProvider.h"
#include "GenericTrackFitter.h"
#include "ThreeDModel.h"

// ROOT includes
//#include "TH1F.h"
//#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"

// Namespaces
using namespace std; 

namespace depfet {
 
/** Constructor 
 */
KalmanAlignmentAlgorithm2::KalmanAlignmentAlgorithm2() {}


/** Jacobian matrix, derivatives of measurement equation f=(fu,fv) to 
 *  the six alignment parameters ar=(dx, dy, dz, dalpha, dbeta, dgamma).
 *
 *  The Jacobian matrix is evaluated for track parematers p0=(tu, tv, u, v)
 *  and for given sensor position and rotation.   
 */
SensorAlignmentJacobian KalmanAlignmentAlgorithm2::Jacobian_Alignment(const TrackState & p0, const Matrix3d & Rot ) const
{
  
  // First, we build the derivatives wrt. the alignment 
  // parameters aq=(du, dv, dw, dalpha, dbeta, dgamma). 
  //  
  // Reference: see V. Karimäki's HIP paper
     
  // Local track parameters  
  double tu = p0[0];
  double tv = p0[1];
  double u  = p0[2];
  double v  = p0[3]; 
  
  // The jacobian matrix
  SensorAlignmentJacobian Jaq;
  
  Jaq(0,0) = -1;      // dfu / ddu
  Jaq(1,0) = 0;       // dfv / ddu
  Jaq(0,1) = 0;       // dfu / ddv
  Jaq(1,1) = -1;      // dfv / ddv
  Jaq(0,2) = +tu;     // dfu / ddw
  Jaq(1,2) = +tv;     // dfv / ddw
  Jaq(0,3) = -v*tu;   // dfu / ddalpha
  Jaq(1,3) = -v*tv;   // dfv / ddalpha
  Jaq(0,4) = u*tu;    // dfu / ddbeta
  Jaq(1,4) = u*tv;    // dfv / ddbeta
  Jaq(0,5) = -v;      // dfu / ddgamma
  Jaq(1,5) = u;       // dfv / ddgamma
   
  //cout << "Jaq " << Jaq << endl; 
  
  // We can apply the chain rule to obtain the Jacobian
  // matrix wrt. ar=(dx, dy, dz, dalpha, dbeta, dgamma).
  
  // Compute A=daq/dar
  Eigen::Matrix<double,6,6> A;
  A = Eigen::Matrix<double,6,6>::Identity();
  A.block<3,3>(0,0) = Rot;
  
  //cout << "A " << A << endl; 
   
  // Apply chain rule 
  auto Ja=Jaq*A; 
  
  //cout << "Ja " << Ja << endl; 
  
  return Ja;
}


/** Align detector with align constants. Returns error flag. 
 */
bool KalmanAlignmentAlgorithm2::AlignDetector(TBDetector& detector, AlignableDet& alignconst)
{
  int nSensors = detector.GetNSensors();  
  
  // Align all sub detectors   
  for (int ipl = 0; ipl < nSensors; ipl++) {
         
    // Load current pixel module    
    Det & adet = detector.GetDet(ipl);
         
    // Read alignment corrections for sensor ipl 
    SensorAlignmentParameters alignPars = alignconst.GetAlignState(ipl);
    double dx = alignPars(0);  
    double dy = alignPars(1);      
    double dz = alignPars(2);     
    double dalpha = alignPars(3); 
    double dbeta  = alignPars(4); 
    double dgamma = alignPars(5);
     
    // Compute a 'delta' frame from corrections 
    ReferenceFrame deltaFrame = ReferenceFrame::create_karimaki_delta(dx,dy,dz,dalpha,dbeta,dgamma); 
    
    // Merge nominal frame and delta frame 
    ReferenceFrame alignFrame = ReferenceFrame::combine_karimaki(adet.GetNominal(), deltaFrame); 
     
    // Update nominal sensor reference frame
    adet.SetNominalFrame(alignFrame); 
    
  }
  
  return false;  
}


/** Performs alignment fit. Returns alignment results. 
 */
AlignableDet KalmanAlignmentAlgorithm2::Fit(TBDetector& detector, TFile * AlignmentData, AlignableDet initialAlignState, 
                   int maxTracks, int annealingTracks, double annealingFactor, double  pValueCut, double deviationCut, bool useBC, int logLevel)
{
 
  AlignmentData->cd("");   
    
   
  // Get information what needs to be aligned
  //------------------------------------------
  int nAlignables = detector.GetNSensors();
  int nParameters = 6;  
  
  cout << "Alignment request for " << nAlignables << " alignable sensors " 
       << " with " << nParameters << " parameters." << endl; 
  
  // Create alignment store 
  //-----------------------
  // It keeps all alignment corrections and quality 
  // indicators. 
  AlignableDet AlignStore = initialAlignState;
   
  
  // Create look up tables of floating alignment variables
  //------------------------------------------------------
  // In general, each alignable (pixel sensor) is considered 
  // to have six rigid body alignment parameters, 3 shifts 
  // and three tilts. 
  // For special reasons, some parameters may be kept fixed in 
  // the alignment procedure. All other alignment paramters  
  // are floating. 
  
  // Fixation of alignables 
  std::vector<bool> fixedAlignable(nAlignables,true);
  
  for (int iAlignable=0; iAlignable < nAlignables; iAlignable++){
    int nFixed = 0; 
    for (int iParameter=0; iParameter < nParameters; iParameter++){
      // Parameter is fixed? 
      double sigma2 =  AlignStore.GetAlignCovariance(iAlignable)(iParameter,iParameter);   
      if (sigma2 == 0 ) {
        ++nFixed;
      }     
    }
    // Do annealing if at least one alignment parameter floating 
    if ( nFixed < 6 ) {
      fixedAlignable[iAlignable] = false; 
      cout << "Alignable " << iAlignable << " is not fixed!" << endl;
    } 
  }
  
  // Get alignment tree 
  TTree * t = (TTree *) AlignmentData->Get("AlignTree");
  if (t == 0) {
    cout << "Error reading object AlignTree from file." << endl;
    return AlignStore;
  }
  
  AlignEvent * alignEvent = new AlignEvent;
  t->SetBranchAddress("AlignEvent", & alignEvent);
  
  // Number of events to be processed (-1 or 0 means all)
  if (maxTracks <= 0)
    maxTracks=t->GetEntriesFast();
  
  cout << "Processing " << maxTracks  << " events" << endl;
  
  // Counters
  int nStep = 0;
  int nOutliers1 = 0;
  int nOutliers2 = 0;
  std::vector<int> nUpdates(nAlignables,0);

  
  // Configure track fitter 
  GenericTrackFitter TrackFitter(detector);
  TrackFitter.SetNumIterations(2);     
  TrackFitter.SetUseBeamConstraint(useBC);
     
  // Loop
  //------
  for (int ii = 0; ii < t->GetEntriesFast(); ii++) {
    
    // Check max events
    if (nStep >= maxTracks) {
      break;
    }
     
    t->GetEntry(ii); 
       
    // This is the nominal detector. Previous tracking was 
    // done in this 'badly aligned' detector
    TBDetector aligned_detector = detector;      
        
    // Next, we update the previous geometry to refit the 
    // track in better aligned detector.  
    AlignDetector(aligned_detector, AlignStore);
    
    // Build TBTrack in aligned detector, i.e. after applying 
    // alignment corrections 
    //--------------------------------
    
    KalmanAlignmentInputProvider kaip;       
    TBTrack trk = kaip.MakeTBTrack( *alignEvent, aligned_detector );  
    
    

    // Deterministic annealing (geometric cooling scheme)
    // must be applied before the track is refitted!!
    // 
    // Annealing multiplies the measurement covariance matrix 
    // of track hits by a scaling factor alpha. The blown 
    // up covariance accounts for the sensor misalignment in 
    // track fitting and alignment.  
    //----------------------------------------------------
     
    for (int ipl= 0; ipl< nAlignables; ++ipl) {  
      
      if (nUpdates[ipl] < annealingTracks) {  
           
        // Get data for alignable 
        //------------------------
        TBTrackElement& TE = trk.GetTE(ipl);  
        
        // Skip sensor w/o measurment
        if ( !TE.HasHit() ) continue;  
        
        // Do not touch fixed alignables 
        if (fixedAlignable[ipl]) continue;  
        
        // This is the actual annealing
        double alpha = TMath::Power(annealingFactor, (annealingTracks-nUpdates[ipl])/(static_cast<float>(annealingTracks)));
        TE.GetHit().GetCov() *= alpha;
        
      } 
    }
    
    
    // ReFit track with (current alignment +  annealing) 
    //------------------------------------------------- 
    bool trkerr = TrackFitter.Fit(trk);
    if ( trkerr ) {
      if (logLevel > 2) { 
        cout << "Fit failed. Skipping track!" << endl;
      }
      continue;
    }  
     
    // Data quality checks I 
    //---------------------
    double prob = TMath::Prob(trk.GetChiSqu(), trk.GetNDF());
    if (prob < pValueCut) { 
      if (logLevel > 2) { 
        cout << "Skipping event because of very low chi2 probability." << endl;
      }
      
      nOutliers1++;
      continue;
    } 
        
    // Kalman update for hit alignables 
    // -------------------------------- 
    // Loop over all alignable sensors
    for (int ipl= 0; ipl< nAlignables; ++ipl) {  
         
      // Get data for alignable 
      //------------------------
      TBTrackElement& TE = trk.GetTE(ipl);  
         
      // Skip sensor w/o measurment
      if ( !TE.HasHit() ) continue;  
      
      // Check module crossing is fitted 
      if ( !TE.IsCrossed() ) continue;
      
      // Do not touch fixed alignables 
      if (fixedAlignable[ipl]) continue; 
      
      // Debug output
      //--------------
      if (logLevel > 2) { 
        cout << "Align update sensor: " << ipl << endl;
      }
      
      // The notation follows Rudi Fruehwirth's paper 
      // "Estimation of Alignment Parameters, using 
      // Kalman Filter with annealing", CMS Note 2002/008.   
      
      // Remember that p0,C0 are track estimates with the 
      // current (global) detector alignment 
       
      auto p0 = TE.GetState().GetPars();
      auto C0 = TE.GetState().GetCov(); 
      auto m = TE.GetHit().GetCoord();     
      auto V = TE.GetHit().GetCov();   
         
      // Copy local alignment state/cov  
      //-------------------------------
      auto a0 = AlignStore.GetAlignState(ipl);
      auto E0 = AlignStore.GetAlignCovariance(ipl); 
      
      // Jacobian matrix, derivatives of measurement equation 
      // m=f(p,a) to track parameters p at (a0,p0) 
      auto H = TrackFitter.GetHMatrix();
      
      // Jacobian matrix, derivatives of measurement equation 
      // m=f(p,a) to alignment parameters a at (a0,p0)  
      auto Rot = TE.GetDet().GetNominal().GetRotation();
      auto D = Jacobian_Alignment(p0, Rot); 
       
      // Weigth matrix of measurment 
      bool invertible = true;
      Matrix2d W = (V + H*C0*H.transpose() + D*E0*D.transpose() ).inverse();
      if (!invertible) {
        if (logLevel > 2) { 
          cout << "ERR: Matrix inversion failed. Skipping alignable!" << endl;
        } 
        continue;
      }	
      
      // Gain matrix K for alignment 
      auto K = E0 * D.transpose() * W;
      // Track prediction for hit coord   
      auto f = H*p0;     
              
      // Update for alignment state + cov  
      auto a1 = a0 + K * (m - f) ;
      auto E1 = E0 - E0*D.transpose()*W*D*E0.transpose();   
          
      // Outlier rejection II
      //----------------------
      bool filterEvent = false;
      if (deviationCut > 0.0){
        for (int ipar=0; ipar < nParameters; ipar++){
          if (TMath::Abs(a1[ipar]-a0[ipar]) > deviationCut * TMath::Sqrt(E0(ipar,ipar))) {
            cout << "Deviation cut for parameter "<< ipar << ": " << TMath::Abs(a1[ipar] - a0[ipar] ) 
                 << " while std. deviation is " << TMath::Sqrt(E0(ipar,ipar)) << endl;
            filterEvent = true;
          }
        }
      }
      
      // Set update for state/cov  
      //-------------------------
      if (!filterEvent){
        AlignStore.SetAlignState(ipl,a1);
        AlignStore.SetAlignCovariance(ipl,E1);    
        ++nUpdates[ipl];
      } else { 
        cout << "Outlier rejection stage 2: Skipping event" << endl;
        nOutliers2++;
      }
         
    } // END UPDATE LOOP   
    
    // Increase alignment step counter
    //---------------------------------
    nStep++;
    // Occasionally give status report
    if (nStep % 10000 == 0){
      TDatime theTime = TDatime();
      cout << theTime.GetHour() << ":" << theTime.GetMinute() << "'" << theTime.GetSecond()
                                << ": " << nStep << " tracks finished ..." << endl;
    }
       
  } // END MAIN LOOP
    
  cout << nOutliers1 << " events rejected by outlier rejection stage 1" << endl;
  cout << nStep << " tracks processed" << endl;
  if (deviationCut > 0.0){
    cout << nOutliers2 << " events ignored by outlier rejection stage 2" << endl;
  } 
  
  return AlignStore; 
}



} // Namespace


