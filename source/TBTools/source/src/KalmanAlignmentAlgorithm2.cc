// KalmanAlignmentAlgorithm2 implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// C++ includes 
#include <iostream>

// DEPFETTrackTool includes
#include "KalmanAlignmentAlgorithm2.h"
#include "AlignEvent.h"
#include "KalmanAlignmentInputProvider.h"
#include "GenericTrackFitter.h"
#include "ThreeDModel.h"

// ROOT includes
#include "TEnv.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"

// Namespaces
using namespace CLHEP;
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
HepMatrix KalmanAlignmentAlgorithm2::Jacobian_Alignment(const HepMatrix & p0, const HepMatrix & Rot,const HepVector & Pos ) const
{
  
  // First, we build the derivatives wrt. the alignment 
  // parameters aq=(du, dv, dw, dalpha, dbeta, dgamma). 
  //  
  // Reference: see V. Karimäki's HIP paper
     
  // Local track parameters  
  double tu = p0[0][0];
  double tv = p0[1][0];
  double u  = p0[2][0];
  double v  = p0[3][0]; 
  
  // The jacobian matrix
  HepMatrix Jaq(2, 6);
  
  Jaq[0][0] = -1;      // dfu / ddu
  Jaq[1][0] = 0;       // dfv / ddu
  Jaq[0][1] = 0;       // dfu / ddv
  Jaq[1][1] = -1;      // dfv / ddv
  Jaq[0][2] = +tu;     // dfu / ddw
  Jaq[1][2] = +tv;     // dfv / ddw
  Jaq[0][3] = -v*tu;   // dfu / ddalpha
  Jaq[1][3] = -v*tv;   // dfv / ddalpha
  Jaq[0][4] = u*tu;    // dfu / ddbeta
  Jaq[1][4] = u*tv;    // dfv / ddbeta
  Jaq[0][5] = -v;      // dfu / ddgamma
  Jaq[1][5] = u;       // dfv / ddgamma
   
  //cout << "Jaq " << Jaq << endl; 
  
  // We can apply the chain rule to obtain the Jacobian
  // matrix wrt. ar=(dx, dy, dz, dalpha, dbeta, dgamma).
  
  // Compute A=daq/dar
  HepMatrix A(6,6,1);  
  A.sub(1,1,Rot); 
  
  //cout << "A " << A << endl; 
   
  // Apply chain rule 
  HepMatrix Ja=Jaq*A; 
  
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
    double dx = alignconst.alignmentParameters[ipl*6+0]; 
    double dy = alignconst.alignmentParameters[ipl*6+1];      
    double dz = alignconst.alignmentParameters[ipl*6+2];    
    double dalpha = alignconst.alignmentParameters[ipl*6+3];
    double dbeta  = alignconst.alignmentParameters[ipl*6+4];
    double dgamma = alignconst.alignmentParameters[ipl*6+5];
     
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
AlignableDet KalmanAlignmentAlgorithm2::Fit(TBDetector& detector, TFile * AlignmentData, string ConfigFile)
{
 
  AlignmentData->cd("");   
   
  // Read config file
  //------------------
  TEnv mEnv(ConfigFile.c_str());
  
  // Number of tracks
  int nMax = mEnv.GetValue("nMaxTracks", 0);
  
  // Set verbosity 
  int LogLevel = mEnv.GetValue("LogLevel", 2);  
  
  // Deterministic annealing
  double AnnealingFactor = mEnv.GetValue("annealingFactor", 1000000.);
  int    AnnealingEvents = mEnv.GetValue("annealingEvents", 0);     

  // Outlier rejection
  double probCut = mEnv.GetValue("probabilityCut", 0.5);
  double deviationCut = mEnv.GetValue("deviationCut", 1.0);
    
  // Use beam model to constrain track fitting
  int useBC = mEnv.GetValue("useBeamModel", 0);
  double beam_divx = mEnv.GetValue("SlopeRMSX", 0.001);
  double beam_divy = mEnv.GetValue("SlopeRMSY", 0.001);
  double beam_sizex = mEnv.GetValue("SpotSizeX", 8.0);
  double beam_sizey = mEnv.GetValue("SpotSizeY", 7.0);
  double beam_corrx = mEnv.GetValue("CorrelationX", 0.0);
  double beam_corry = mEnv.GetValue("CorrelationY", 0.0);
  
  // Print info about alignment options
  //------------------------------------
  if (probCut > 0.0) { cout << "Track probability cut: " << probCut << endl; }
  if (deviationCut > 0.0) { cout << "Deviation cut: " << deviationCut << endl; }
  if (AnnealingEvents != 0) { cout << "Deterministic annealing until event " << AnnealingEvents << endl; }
  if (AnnealingEvents != 0) { cout << "Initial annealing factor is " << AnnealingFactor << endl; }
  if (useBC) { 
    cout << "Using beam constraint: " << endl; 
    cout << " beam spot sizes [mm]: " << beam_sizex << ", " << beam_sizey << endl; 
    cout << " beam slope rms [rad]: " << beam_divx  << ", " << beam_divy  << endl; 
    cout << " beam correlations: " << beam_corrx  << ", " << beam_corry  << endl; 
  }

  
  
  // Print info about alignment options
  //------------------------------------
  if (probCut > 0.0) { cout << "Track p-value cut: " << probCut << endl; }
   
  
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
  AlignableDet AlignStore(nAlignables,nParameters);
   
  // Initialize parameters and covariance
  //--------------------------------------
  for (int iAlignable=0; iAlignable < nAlignables; iAlignable++){
    for (int iParameter=0; iParameter < nParameters; iParameter++){
      double error =  mEnv.GetValue(Form("errorParameter%dofAlignable%d", iParameter, iAlignable ), 0.0); 
      
      cout << "Set parameter " << iParameter << " of alignable " << iAlignable  << " to 0"
           << " +- " << error << endl;
              
      AlignStore.alignmentParameters[iAlignable*nParameters + iParameter] = 0;
      AlignStore.alignmentCovariance[iAlignable*nParameters + iParameter][iAlignable*nParameters + iParameter] = error*error;
    }
  }
  
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
      double error =  mEnv.GetValue(Form("errorParameter%dofAlignable%d", iParameter, iAlignable ), 0.0); 
      if (error == 0 ) {
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
  if (nMax <= 0)
    nMax=t->GetEntriesFast();
  
  cout << "Processing " << nMax  << " events" << endl;
  
  
  // Counters
  int nStep = 0;
  int nOutliers1 = 0;
  int nOutliers2 = 0;
  std::vector<int> nUpdates(nAlignables,0);

  
  // Configure track fitter 
  GenericTrackFitter TrackFitter(detector);
  TrackFitter.SetNumIterations(2);     
  TrackFitter.SetUseBeamConstraint(useBC);
  TrackFitter.SetBeamDivergenceX(beam_divx);
  TrackFitter.SetBeamDivergenceY(beam_divy);
  TrackFitter.SetBeamSizeX(beam_sizex);
  TrackFitter.SetBeamSizeY(beam_sizey);
  TrackFitter.SetBeamCorrelationX(beam_corrx);
  TrackFitter.SetBeamCorrelationY(beam_corry);
  
     
  // Loop
  //------
  for (int ii = 0; ii < t->GetEntriesFast(); ii++) {
    
    // Check max events
    if (nStep >= nMax) {
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
      
      if (nUpdates[ipl] < AnnealingEvents) {  
           
        // Get data for alignable 
        //------------------------
        TBTrackElement& TE = trk.GetTE(ipl);  
        
        // Skip sensor w/o measurment
        if ( !TE.HasHit() ) continue;  
        
        // Do not touch fixed alignables 
        if (fixedAlignable[ipl]) continue;  
        
        // This is the actual annealing
        double alpha = TMath::Power(AnnealingFactor, (AnnealingEvents-nUpdates[ipl])/(static_cast<float>(AnnealingEvents)));
        TE.GetHit().GetCov() *= alpha;
        
      } 
    }
    
    
    // ReFit track with (current alignment +  annealing) 
    //------------------------------------------------- 
    bool trkerr = TrackFitter.Fit(trk);
    if ( trkerr ) {
      if (LogLevel > 2) { 
        cout << "Fit failed. Skipping track!" << endl;
      }
      continue;
    }  
     
    // Data quality checks I 
    //---------------------
    double prob = TMath::Prob(trk.GetChiSqu(), trk.GetNDF());
    if (prob < probCut) { 
      if (LogLevel > 2) { 
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
      if (LogLevel > 2) { 
        cout << "Align update sensor: " << ipl << endl;
      }
      
      // The notation follows Rudi Fruehwirth's paper 
      // "Estimation of Alignment Parameters, using 
      // Kalman Filter with annealing", CMS Note 2002/008.   
      
      // Remember that p0,C0 are track estimates with the 
      // current (global) detector alignment 
       
      HepMatrix p0 = TE.GetState().GetPars();
      HepSymMatrix C0 = TE.GetState().GetCov(); 
      HepMatrix m = TE.GetHit().GetCoord();     
      HepSymMatrix V = TE.GetHit().GetCov();   
         
      // Copy local alignment state/cov  
      //-------------------------------
      HepVector a0 = AlignStore.GetAlignState(ipl);
      HepSymMatrix E0 = AlignStore.GetAlignCovariance(ipl); 
      
      // Jacobian matrix, derivatives of measurement equation 
      // m=f(p,a) to track parameters p at (a0,p0) 
      HepMatrix H = TrackFitter.GetHMatrix();
      
      // Jacobian matrix, derivatives of measurement equation 
      // m=f(p,a) to alignment parameters a at (a0,p0) 
      HepVector Pos = TE.GetDet().GetNominal().GetPosition(); 
      HepMatrix Rot = TE.GetDet().GetNominal().GetRotation();
      HepMatrix D = Jacobian_Alignment(p0, Rot, Pos); 
       
      // Weigth matrix of measurment 
      int ierr; 
      HepSymMatrix W = (V + C0.similarity(H) + E0.similarity(D) ).inverse(ierr);
      if (ierr != 0) {
        if (LogLevel > 2) { 
          cout << "ERR: Matrix inversion failed. Skipping alignable!" << endl;
        } 
        continue;
      }	
      
      // Gain matrix K for alignment 
      HepMatrix K = E0 * D.T() * W;
      // Track prediction for hit coord   
      HepMatrix f = H*p0;     
              
      // Update for alignment state + cov  
      HepVector a1 = a0 + K * (m - f) ;
      HepSymMatrix E1 = E0 - (W.similarityT(D)).similarity(E0); 
       
      // Outlier rejection II
      //----------------------
      bool filterEvent = false;
      if (deviationCut > 0.0){
        for (int ipar=0; ipar < nParameters; ipar++){
          if (TMath::Abs(a1[ipar]-a0[ipar]) > deviationCut * TMath::Sqrt(E0[ipar][ipar])) {
            cout << "Deviation cut for parameter "<< ipar << ": " << TMath::Abs(a1[ipar] - a0[ipar] ) 
                 << " while std. deviation is " << TMath::Sqrt(E0[ipar][ipar]) << endl;
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


