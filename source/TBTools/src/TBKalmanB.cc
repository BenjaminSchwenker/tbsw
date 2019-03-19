// TBKalmanB implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Marlin includes 
#include <marlin/Global.h>
#include <streamlog/streamlog.h>
  
// TBTools includes
#include "MaterialEffect.h"
#include "TBKalmanB.h"
#include "HelixTrackModel.h"
#include "StraightLineTrackModel.h"
#include "ThreeDModel.h" 

// C++ STL includes 
#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/LU>

// Namespaces
using namespace marlin;
using namespace std; 

namespace depfet {



/** Constructor 
 */
TBKalmanB::TBKalmanB(TBDetector& detector)
{ 
  // Number of iterations for double filter
  NumIt = 1;       
    
  // Project hit coord out of track state
  H = StateHitProjector::Zero();
  H(0,2) = 1.;
  H(1,3) = 1.;

  // Set the initial track state vector
  x0_init = TrackState::Zero();
    
  // Set the initial track state covariance matrix
  C0_init = TrackStateCovariance::Zero();
    
  for (int i = 0; i < 2; i++) {
    // There is a tradeoff between choosing large values (which makes the
    // tracking numerically unstable) and choosing small values (which biases
    // the track from the seed). I found 1E4 working for this model, but this
    // might need to get tuned.
    C0_init(i,i) = 1E-2;
    C0_init(i+2,i+2) = 1E2; 
  } 
  C0_init(4,4) = 1;  
  
  // Electron mass, in GeV 
  mass = 0.000511;

  // Electron charge in, e
  charge = -1;
  
  // TODO: In order to use the beam constraint, the inital covariance 
  // matrix neeed to be computed. this code is missing rigth now. 
  // Currently, using the useBC flag will be ignored.
  
  // By default, beam constraint is not used 
  m_useBC = false; 
  m_sizex = 0;       // spot size, mm
  m_sizey = 0;       // spot size, mm
  m_divx = 0;     // divergence, rad
  m_divy = 0;     // divergence, rad
  m_corrx = 0;      // beam correlation coeff  
  m_corry = 0;      // beam correlation coeff
  
  Vector3d field;
  field<< detector.GetBx(),detector.GetBy(),detector.GetBz();
  
  if ( field.norm() == 0 ) {
    TrackModel = new StraightLineTrackModel();     
  } else {
    TrackModel = new HelixTrackModel(field); 
  }
}

/** Destructor 
 */
TBKalmanB::~TBKalmanB()
{ 
  delete TrackModel;  
}

/** Extrapolates the track seed to all planes. Returns error flag. 
*/
bool TBKalmanB::ExtrapolateSeed(TBTrack& trk) 
{
  
  // Get particle hypothesis
  mass = trk.GetMass();
  charge = trk.GetCharge();
   
  // Linearize fit around a reference track state  
  TrackState RTrkState = trk.GetReferenceState().GetPars();      
  ReferenceFrame RTrkSurf = trk.GetTE( trk.GetReferenceState().GetPlane() ).GetDet().GetNominal();
  
  // Compute intersections of the reference trajectory
  // with sub detectors. 
  
  // Loop over all track elements 
  int nTE = trk.GetNumTEs();  
    
  for(int iTE=0; iTE<nTE; ++iTE) {
    
    // Get track element 
    TBTrackElement& TE = trk.GetTE( iTE ); 
      
    // Next surface along beam line  
    const ReferenceFrame& Next_Surf = TE.GetDet().GetNominal();
           
    // Extrapoate reference state 
    bool error = false;
    TE.GetState().GetPars() = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);  
      
    if (error) {
      TE.SetCrossed(false);
    } else {
      TE.SetCrossed(true); 
    }
    
  } 

  return false;
}


/** Performs track fitting. Returns error flag.
 */
bool TBKalmanB::Fit(TBTrack& trk, int smoothPlane)
{
   
  // Particle hypothesis  
  mass = trk.GetMass();
  charge = trk.GetCharge();
  
  // Number of sensor planes
  int nTE = trk.GetNumTEs();
   
  // Number of planes crossed by seed trajectory
  int nCrossed = 0;
  
  // List of crossed track elements. Predefined here, so it can be reused.
  vector<int> CrossedTEs;
  CrossedTEs.reserve(nTE);
  // Reference trajectory.  Predefined here, so it can be reused.
  vector<TrackState> RefStateVec;
  RefStateVec.reserve(nTE);
  // List of predicted estimators from forward pass.  Predefined here, so it can be reused.
  FilterStateVec FFilterData;
  FFilterData.reserve(nTE);
  // List of predicted estimators from backward pass.  Predefined here, so it can be reused.
  FilterStateVec BFilterData;  
  BFilterData.reserve(nTE);
  
  // Improve the seed parameters used for linearization of the fit. 
  // 
  // We can improve the seed by running a backward
  // filter from the last to the firts plane in the beam. 
  // A new seed is taken to be the predicted state at 
  // the first plane. 
  
  for(int iter=0; iter<NumIt; iter++) { 
    
    // TODO: Extrapolation of the reference state into the telescope is a repeating topic. One could deal with it 
    // by reusing the ExtrapolateSeed() member function. But one needs to circumvent the usage of the CrossedTEs
    // vector in the fitting code. 
    
    // Get paramters of reference trajectory   
    TrackState RTrkState = trk.GetReferenceState().GetPars();      
    ReferenceFrame RTrkSurf = trk.GetTE( trk.GetReferenceState().GetPlane() ).GetDet().GetNominal(); 
     
    // Clear list of crossed track elements
    CrossedTEs.clear();
    
    // Clear reference trajectory
    RefStateVec.clear();
    
    // Loop over all track elements (sensor planes)
    for(int iTE=0; iTE<nTE; ++iTE) {
      
      // Next surface along beam line  
      const ReferenceFrame& Next_Surf = trk.GetTE(iTE).GetDet().GetNominal();
           
      // Extrapoate reference state 
      bool error = false;
      TrackState Next_State = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);  
      
      if (error) {
        trk.GetTE(iTE).SetCrossed(false);
      } else {
        trk.GetTE(iTE).SetCrossed(true);
        // Remember plane number 
        CrossedTEs.push_back(iTE);
        // Remember reference state 
        RefStateVec.push_back(Next_State);      
      }
    } 
    
    // Number of crossed track elements 
    nCrossed = (int) CrossedTEs.size(); 
    
    // The reference trajectory should at least 
    // intersect with one detector.  
    
    if ( nCrossed == 0 ) {
      trk.SetChiSqu(-1);
      SetNdof(trk);
      return true;
    }   
    
    // Init backward filter 
    // The backward pass is initialized using 
    // the reference track state at the last
    // sensor. 
    BFilterData.clear();
    BFilterData[nCrossed-1].Pr_x = x0_init;   
    BFilterData[nCrossed-1].Pr_C = C0_init;  
    
    // Backward FILTER  -----------------------------------------------------
    double chisqu_back = FilterPass(trk, CrossedTEs, RefStateVec, -1, BFilterData);  
      
    // A very basic consistency test 
    if ( std::isnan(chisqu_back) ||  chisqu_back < 0 )  {
      trk.SetChiSqu(-1);
      SetNdof(trk);
      return true;
    } 
    
    // Set new seed parameters and new linearization point  
    trk.GetReferenceState().Pars =  RefStateVec[0] +  BFilterData[0].Pr_x;   
    trk.GetReferenceState().SetPlane(CrossedTEs[0]);  
    trk.SetMomentum(std::abs(charge/trk.GetReferenceState().Pars[4]));
  } // end iterations
   

  // Recompute the seed trajectory 
  //
  // Get the reference parameters and extrapolate them into 
  // the telescope. Same as before. 
  
  // Get paramters of reference trajectory   
  TrackState RTrkState = trk.GetReferenceState().GetPars();      
  ReferenceFrame RTrkSurf = trk.GetTE( trk.GetReferenceState().GetPlane() ).GetDet().GetNominal();
  

// TODO: This code is currently brocken and should not be used.     
//  if ( m_useBC ) {
    // The idea is to use the beam constraint as a reference
    // trajectory 
    
//    double mom = trk.GetMomentum();    
//    RTrkSurf = trk.GetTE(0).GetDet().GetNominal();
//    RTrkState = ComputeBeamConstraint(RTrkSurf, mom, charge);
//  }  
    
  // Clear list of crossed track elements
  CrossedTEs.clear();
    
  // Clear reference trajectory
  RefStateVec.clear();
    
  // Compute intersections of the reference trajectory
  // with sub detectors. 
  
  // Loop over all sensor planes
  for(int iTE=0; iTE<nTE; ++iTE) {
      
    // Next surface along beam line  
    const ReferenceFrame& Next_Surf = trk.GetTE(iTE).GetDet().GetNominal();
           
    // Extrapoate reference state 
    bool error = false;
    TrackState Next_State = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);  
      
    if (error) {
      trk.GetTE(iTE).SetCrossed(false);
    } else {
      trk.GetTE(iTE).SetCrossed(true);
      // Remember plane number 
      CrossedTEs.push_back(iTE);
      // Remember reference state 
      RefStateVec.push_back(Next_State);      
    }
  } 
    
  // Number of crossed track elements  
  nCrossed = (int) CrossedTEs.size(); 
    
  // The reference trajectory should at least 
  // intersect with one detector.  
    
  if ( nCrossed == 0 ) {
    trk.SetChiSqu(-1);
    SetNdof(trk);
    return true;
  }   
     
  // Init forward filter 
  // The forward pass is initialized using 
  // the reference track state at the first
  // sensor. 
  FFilterData.clear();  
  FFilterData[0].Pr_x = x0_init;
  FFilterData[0].Pr_C = C0_init;
    
  // Run forward FILTER -----------------------------------------------------
  double chisqu = FilterPass(trk, CrossedTEs, RefStateVec, +1, FFilterData); 
  
  // A very basic consistency test 
  if ( std::isnan(chisqu) ||  chisqu < 0 )  {
    trk.SetChiSqu(-1);
    SetNdof(trk);
    return true;
  }
  
  // Init backward filter 
  // The backward pass is initialized using 
  // the reference track state at the last
  // sensor. 
  BFilterData.clear();
  BFilterData[nCrossed-1].Pr_x = x0_init;   
  BFilterData[nCrossed-1].Pr_C = C0_init;  
    
  // Run backward FILTER  -----------------------------------------------------
  double chisqu_back = FilterPass(trk, CrossedTEs, RefStateVec, -1, BFilterData);  
      
  // A very basic consistency test 
  if ( std::isnan(chisqu_back) ||  chisqu_back < 0 )  {
    trk.SetChiSqu(-1);
    SetNdof(trk);
    return true;
  } 

  // Smoothing 
  // 
  // Finally, we carry out the smoothing of forwars backward 
  // filter estimates. This is done by combining the forward 
  // and backward estimates at all crossed track elements. 
       
  for(int is=0; is<nCrossed; ++is) {
       
    // Get plane number 
    int iTE = CrossedTEs[is];  
    
    // Smoothing is expensive and should only be done for 
    // all planes when requested by the user. 
    if ( (smoothPlane < 0) || (iTE == smoothPlane)  ) { 
      
      // Get track element 
      TBTrackElement& TE = trk.GetTE( iTE ); 
      
      // Compute unbiased smoothed estimate 
      // of local track state. 
          
      TrackState & xs =TE.GetState().GetPars(); 
      TrackStateCovariance & Cs = TE.GetState().GetCov();
        
      // Smoothing 
      if (is == 0) {
        xs = BFilterData[is].Pr_x; 
        Cs = BFilterData[is].Pr_C;  
      } else if (is == nCrossed-1) {
        xs = FFilterData[is].Pr_x;    
        Cs = FFilterData[is].Pr_C;
      } else {
        const auto& xb = BFilterData[is].Pr_x;
        const auto& Cb = BFilterData[is].Pr_C;
        const auto& xf = FFilterData[is].Pr_x;
        const auto& Cf = FFilterData[is].Pr_C;
        
        bool error = GetSmoothedData(xb, Cb, xf, Cf, xs, Cs);
        if ( error ) {
          trk.SetChiSqu(-1);
          SetNdof(trk); 
          return true; 
        }
      }
      
      // Linearization: This adds the reference parameters         
      xs += RefStateVec[is];
      
      // Compute local ChiSqu using the smoothed state  
      if ( TE.HasHit() ) {     
        TE.SetChiSqu(GetPredictedChi2(xs, Cs, TE.GetHit()));
      }
    } 
  }  
  
  // Everything is ok
  trk.SetChiSqu(chisqu);
  SetNdof(trk);  
  return false;
}

/** Filters a new hit 
 */
double TBKalmanB::FilterHit(const TBHit& hit, const TrackState& xref, TrackState& x0, TrackStateCovariance& C0) 
{ 
  // Here, the measurment update takes place  
  double predchi2 = 0; 
    
  // Measured hit coordinates, 2x1 matrix 
  const Vector2d& m = hit.GetCoord();
        
  // Covariance for hit coordinates, 2x2 matrix 
  const Matrix2d& V= hit.GetCov(); //CHECK WHY Problematic

  bool invertible = true;
  Matrix2d W = Matrix2d::Zero();// = (V + H*C0*H.transpose()).inverse();
  (V+H*C0*H.transpose()).computeInverseWithCheck(W,invertible);   // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
  if (!invertible) {
    streamlog_out(MESSAGE3) << "Hit filtering: Matrix inversion failed! Returns chi2 < 0!"
                         << std::endl;
    return -1;
  }	
     
  // This is the predicted residual 
  Vector2d r = m - H*(x0 + xref); 
      
  // This is the predicted chi2 
  predchi2 = (r.transpose()*W*r)[0];
      
  // Kalman gain matrix K 
  auto K = C0 * H.transpose() * W; 
       
  // This is the filtered state
  x0 += K * r;
  C0 -= ( C0*H.transpose()*W*H*C0.transpose() ).eval() ; 

        
  return predchi2;
}


/** Propagte state from trackelement te to next track element nte 
 */
int TBKalmanB::PropagateState(TBTrackElement& te, TBTrackElement& nte, TrackState& xref, TrackState& nxref, TrackState& x0, TrackStateCovariance& C0)
{ 
       
  // Reference frame for next track element          
  const ReferenceFrame& nSurf = nte.GetDet().GetNominal();

  // Reference frame for current track element          
  const ReferenceFrame& Surf = te.GetDet().GetNominal();
  
  // Direction of propagation (idir > 0 means along beam direction)
  int idir = nte.GetDet().GetPlaneNumber() - te.GetDet().GetPlaneNumber();
   
  // This fitter takes into account scatter in air
  // gaps between sensors. Therefor, we add a virtual 
  // air surface between the detectors which scatters 
  // the track with a material budget equal to the 
  // extrapolation step length between the sensors.
        
  // Get signed flight length in air between detectors 
  double length = TrackModel->GetSignedStepLength(xref, Surf, nSurf); 
   

  // Extraploate half step along straight line 
  TrackState xref_air = xref;
  ReferenceFrame Surf_air = Surf; 
  TrackModel->Extrapolate(xref_air, Surf_air, length/2);
    
  // To start ierr is set to 0 (= OK)
  int ierr = 0; 

  if ( idir > 0 ) {
          
    // MAP estimate [x0,C0] from te to air surface
    // ---------------------------------------
    double l0 = te.GetDet().GetTrackLength(xref[2], xref[3], xref[0], xref[1]);
    double X0 = te.GetDet().GetRadLength(xref[2],xref[3]);
    double theta2_det = materialeffect::GetScatterTheta2(xref, l0, X0, mass, charge );   
    ierr = MAP_FORWARD( theta2_det, xref, Surf, Surf_air, x0, C0 );
    if (ierr != 0) {
      streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
      return -1;
    }	
         
    // MAP estimate [x0,C0] from air to next te 
    // ---------------------------------------
    double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge);   
    ierr = MAP_FORWARD( theta2_air, xref_air, Surf_air, nSurf, x0, C0 );
    if (ierr != 0) {
      streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
      return -1;
    }	
          
  } else {
        
    // MAP estimate [x0,C0] from old det to air surface
    // ---------------------------------------
    double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge) ;   // Backward form 
    ierr = MAP_BACKWARD( theta2_air, xref, Surf, xref_air, Surf_air, x0, C0 );
    if (ierr != 0) {
      streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
      return -1;
    }	
        
    // MAP estimate [x0,C0] from air to new det surface
    // ---------------------------------------
    double l0 = nte.GetDet().GetTrackLength(nxref[2], nxref[3], nxref[0], nxref[1]);
    double X0 = nte.GetDet().GetRadLength(nxref[2],nxref[3]);
    double theta2_det = materialeffect::GetScatterTheta2(xref, l0, X0, mass, charge);   // Backward form    
    ierr = MAP_BACKWARD( theta2_det, xref_air, Surf_air, nxref, nSurf, x0, C0 );          
    if (ierr != 0) {
      streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
      return -1;
    }	
               
  }
      
  return 0;
}


/** Runs filter. Returns fit chi2. 
 */
double TBKalmanB::FilterPass(TBTrack& trk, std::vector<int>& CrossedTEs, std::vector<TrackState>& RefStateVec, int idir, FilterStateVec& Result) 
{ 
    
  // To start ierr is set to 0 (= OK)
  int ierr = 0; 
  // Final chisqu of track fit
  double chisqu = 0;
       
  std::vector<TBTrackElement>& TEVec = trk.GetTEs();  
  
  // Loop over crossed track elements 
   
  int nCrossed = (int) CrossedTEs.size(); 
  int istart, istop;   
  if (idir > 0) {
    istart = 0; 
    istop = nCrossed-1; 
  } else {
    istart = nCrossed-1; 
    istop = 0; 
  } 
  
  for(int is=istart; is!=istop+idir; is+=idir) {
    
    // Get plane number 
    int iTE = CrossedTEs[is];  
    
    // Process next track element 
    TBTrackElement& te = TEVec[iTE];
    
    // Get surface parameters           
    const ReferenceFrame& Surf = te.GetDet().GetNominal();
    
    // Get current reference state xref
    const auto& xref = RefStateVec[is];

    // Copy predicted [x,C]
    TrackState x0 = Result[is].Pr_x;  
    TrackStateCovariance C0 = Result[is].Pr_C;
     
    // Here, the measurment update takes place  
    double predchi2 = 0; 
    
    if ( te.HasHit() ) {
      
      // Measured hit coordinates, 2x1 matrix 
      const Vector2d& m = te.GetHit().GetCoord();
        
      // Covariance for hit coordinates, 2x2 matrix 
      const Matrix2d& V = te.GetHit().GetCov();
            
      // Weigth matrix of measurment 
      bool invertible = true;
      Matrix2d W = Matrix2d::Zero(); // = (V + H*C0*H.transpose()).inverse();
      (V + H*C0*H.transpose()).computeInverseWithCheck(W,invertible);  // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
      if (!invertible) {
        streamlog_out(MESSAGE3) << "Hit filtering: Matrix inversion failed. Quit filter pass!"
                             << std::endl;
        return -1;
      }	
       
      // This is the predicted residual 
      Vector2d r = m - H*(x0 + xref); 
      
      // This is the predicted chi2  
      predchi2 = (r.transpose()*W*r) [0];   
      
      // Kalman gain matrix K 
      auto K = C0 * H.transpose() * W; 
       
      // This is the filtered state
      x0 += K * r;
      C0 -= ( C0*H.transpose()*W*H*C0.transpose() ).eval();       
    }  
    
    // Store results and update chisqu
    chisqu += predchi2;
     
    // Propagate to next surface
    int inext = is+idir;
  
    if (inext!=istop+idir)  {
         
      // Next surface along filter direction 
      TBTrackElement& nte = TEVec[ CrossedTEs[inext] ];
      
      // Parameters for next surface          
      const ReferenceFrame& nSurf = nte.GetDet().GetNominal();

      // Get reference state  
      TrackState & nxref = RefStateVec[inext];  

      // This fitter takes into account scatter in air
      // gaps between sensors. Therefor, we add a virtual 
      // air surface between the detectors which scatters 
      // the track with a material budget equal to the 
      // extrapolation step length between the sensors.
        
      // Extrapolate to air surface between detector planes 
      TrackState xref_air = xref; 
      ReferenceFrame Surf_air = Surf; 
      
      // Get signed flight length in air between detectors 
      double length = TrackModel->GetSignedStepLength(xref, Surf, nSurf); 
      
      // Extraploate half step along straight line 
      TrackModel->Extrapolate(xref_air, Surf_air, length/2);
      
      if ( idir > 0 ) {
          
        // MAP estimate [x0,C0] from old det to air surface
        // ---------------------------------------
        double l0 = te.GetDet().GetTrackLength(xref[2], xref[3], xref[0], xref[1]);
        double X0 = te.GetDet().GetRadLength(xref[2],xref[3]);    
        double theta2_det = materialeffect::GetScatterTheta2(xref, l0, X0, mass, charge);   
        ierr = MAP_FORWARD( theta2_det, xref, Surf, Surf_air, x0, C0 );
        if (ierr != 0) {
          streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
         
        // MAP estimate [x0,C0] from air to new det surface
        // ---------------------------------------
        double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge);   
        ierr = MAP_FORWARD( theta2_air, xref_air, Surf_air, nSurf, x0, C0 );
        if (ierr != 0) {
          streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
        
        
      } else {
        
        // MAP estimate [x0,C0] from old det to air surface
        // ---------------------------------------
        double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge) ;   // Backward form 
        ierr = MAP_BACKWARD( theta2_air, xref, Surf, xref_air, Surf_air, x0, C0 );
        if (ierr != 0) {
          streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
        
        // MAP estimate [x0,C0] from air to new det surface
        // ---------------------------------------
        double l0 = nte.GetDet().GetTrackLength(nxref[2], nxref[3], nxref[0], nxref[1]);
        double X0 = nte.GetDet().GetRadLength(nxref[2],nxref[3]);    
        double theta2_det = materialeffect::GetScatterTheta2(xref, l0, X0, mass, charge );   // Backward form    
        ierr = MAP_BACKWARD( theta2_det, xref_air, Surf_air, nxref, nSurf, x0, C0 );          
        if (ierr != 0) {
          streamlog_out(MESSAGE2) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
         	
         
      }
      
      // Store results        
      Result[inext].Pr_x = x0;  
      Result[inext].Pr_C = C0;  
      
    } // end propagation 
     
           
  }
    
  
  // Successfull return  
  return chisqu;  
}

/** Compute the weighted means of forward and backward filter estimated. In a 
 *  numerically robust way. Returns error flag.
 */
bool TBKalmanB::GetSmoothedData( const TrackState& xb, const TrackStateCovariance& Cb, const TrackState& xf, const TrackStateCovariance& Cf,
               TrackState& Smoothed_State, TrackStateCovariance& Smoothed_Cov){
    // Compute forward weight matrix
    TrackStateCovariance Cb_inv=Cb.llt().solve(TrackStateCovariance::Identity());
    if(!(Cb* (Cb_inv*xb)).isApprox(xb,1e-6)){
        streamlog_out(MESSAGE3)<< "Smoothing 1: Matrix inversion failed. Quit fitting!"
                 << std::endl;
//        streamlog_out(MESSAGE3)<< "Smoothing 1: Original matrix" << std::endl<< Cb << std::endl;
//        streamlog_out(MESSAGE3)<< "Smoothing 1: Inverted matrix" << std::endl<< Cb_inv << std::endl;
//        streamlog_out(MESSAGE3)<< "Smoothing 1: Other inverse function" << std::endl<< Cb.inverse() << std::endl;
//        streamlog_out(MESSAGE3)<< "Smoothing 1: Multiply result" << std::endl<< Cb*Cb_inv << std::endl;
//        streamlog_out(MESSAGE3)<< "Smoothing 1: Vector" << std::endl<< xb.transpose() << std::endl;
//        streamlog_out(MESSAGE3)<< "Smoothing 1: Vector result" << std::endl<< (Cb* (Cb_inv*xb)).transpose() << std::endl;
//        streamlog_out(MESSAGE3)<< "Smoothing 1: Diff " << std::endl<< (xb-Cb* (Cb_inv*xb)).transpose() << std::endl;
        return true;
    }
    // Compute backward weight matrix
    TrackStateCovariance Cf_inv=Cf.llt().solve(TrackStateCovariance::Identity());
    if(!(Cf*(Cf_inv*xf)).isApprox(xf,1e-6)){
        streamlog_out(MESSAGE3)<< "Smoothing 2: Matrix inversion failed. Quit fitting!"
                 << std::endl;
        return true;
    }
    // Weighted (smoothed) covariance
    Smoothed_Cov=(Cf_inv+Cb_inv).llt().solve(TrackStateCovariance::Identity());
    // Weighted state
    Smoothed_State = Smoothed_Cov*(Cf_inv*xf + Cb_inv*xb);
    if(!((Cf_inv+Cb_inv)* Smoothed_State).isApprox(Cf_inv*xf + Cb_inv*xb,1e-6)){
        streamlog_out(MESSAGE3)<< "Smoothing 3: Matrix inversion failed. Quit fitting!"
                 << std::endl;
        return true;
    }
    return false;
}
/* Old implementation*/
//bool TBKalmanB::GetSmoothedData( TrackState& xb, TrackStateCovariance& Cb, TrackState& xf, TrackStateCovariance& Cf,
//                                TrackState& Smoothed_State, TrackStateCovariance& Smoothed_Cov)
//{

//  // Error flag
//  bool invertible = true;

//  // Compute forward weight matrix
//  TrackStateWeight fW = Cf.inverse();
//  //Cf.computeInverseWithCheck(fW,invertible);
//  if (!invertible) {
//    streamlog_out(MESSAGE3) << "Smoothing 1: Matrix inversion failed. Quit fitting!"
//                         << std::endl;
//    return invertible;
//  }

//  // Compute backward weight matrix
//  TrackStateWeight bW = Cb.inverse();
//  //Cb.computeInverseWithCheck(bW,invertible);
//  if (!invertible) {
//    streamlog_out(MESSAGE3) << "Smoothing 2: Matrix inversion failed. Quit fitting!"
//                         << std::endl;
//    return invertible;
//  }

//  // Weighted (smoothed) covariance
//  Smoothed_Cov = (fW + bW).inverse();
//  //(fW + bW).computeInverseWithCheck(Smoothed_Cov,invertible);
//  if (!invertible) {
//    streamlog_out(MESSAGE3) << "Smoothing 3: Matrix inversion failed. Quit fitting!"
//                         << std::endl;
//    return invertible;
//  }

//  // Weighted state
//  Smoothed_State = Smoothed_Cov*(fW*xf + bW*xb);

//  // Successfull return
//  return !invertible;
//}

int TBKalmanB::MAP_FORWARD(  double theta2, 
                        const TrackState& xref, const ReferenceFrame& Surf, const ReferenceFrame& nSurf,
                        TrackState& x0, TrackStateCovariance& C0
                     )
{ 
  
  // The linearized mapping from state at plane k to state at plane k+1:
  // 
  //  x_k+1 = J(k+1,k) * x_k + G_k * w_k                        (*)
  // 
  //  J(k+1,k): Transport matrix from k to k+1   
  //  w_k     : Vector of scatter angles  at k 
  //  x_k     : Vector of track parameters at k
  //  Q_k     : Diag. covariance matrix of w_k 
  //  G_k     : Scatter gain matrix at k    
  // 
  // The transport matrix J(k+1,k) and the scatter gain matrix G_k are computed
  // for the reference track.  
  // 
  // The first term in equation (*) is just a straight line extrapolation, 
  // but the second term includes the influence of random scatterings w_k
  // at plane k. From a physics perspective, x_k is a 'in' state, i.e. 
  // the vector of track parameters before the scattering at plane k happens.
  // 
  // For the state extrapolation in the filter, we can only use the a priori
  // info about scattering angles, i.e <w_k> = 0. For the covariance update, 
  // we can estimate the amount of added uncertainty due to multiple scatter
  // noise.
        
  // Time update of covariance  matrix  
  int ierr = 0; 
  TrackStateJacobian J;      
  ierr = TrackModel->TrackJacobian( xref, Surf, nSurf, J);  
  if (ierr != 0) {
    return -1;
  }	

  // Add scatter noise
  // -----------------------------------
        
  // Local Scatter gain matrix      
  TrackStateGain Gl = TrackModel->GetScatterGain(xref);
            
  // General Scatter gain matrix 
  TrackStateGain G = J*Gl;   
  
  x0 = (J*x0).eval();
  C0.noalias() = (J*C0).eval()*J.transpose() + theta2*G*G.transpose();

  return 0; 
}

int TBKalmanB::MAP_BACKWARD(  double theta2, 
                        const TrackState& xref, const ReferenceFrame& Surf,
                        const TrackState& nxref, const ReferenceFrame& nSurf,
                        TrackState& x0, TrackStateCovariance& C0
                     )
{
   
  // The linearized mapping from state x_k+1 at plane k+1 to state x_k at plane k:
  // 
  // x_k = J(k+1,k)^-1 * x_k+1 - J(k+1,k)^-1 * G_k * w_k   (**) 
  // 
  //  J(k+1,k): Transport matrix from k to k+1   
  //  w_k     : Vector of scatter angles at k 
  //  x_k     : Vector of track parameters at k
  //  Q_k     : Diag. covariance matrix of w_k 
  //  G_k     : Scatter gain matrix at k   
  // 
  // The transport matrix J(k+1,k) and the scatter gain matrix G_k are computed
  // for the reference track.  
  // 
  // The first term in equation (**) is just a straight line extrapolation from 
  // plane k+1 to plane k. This extrapolation uses the inverse transport matrix 
  // J(k+1,k)^-1 = J(k,k+1). The second term includes the influence of random 
  // scatterings w_k at plane k. Note that x_k depends only on the scatterings w_k 
  // at plane k and not(!!) on the scatterings at plane k+1.  From a physics perspective,
  // x_k is a 'in' state, i.e. the vector of track parameters before the scattering
  // at plane k happens.
  // 
  // For the state extrapolation in the filter, we can only use the a priori
  // info about scattering angles, i.e <w_k> = 0. For the covariance update, 
  // we can estimate the amount of added uncertainty due to multiple scatter
  // noise.
  
  // Time update of covariance  matrix 
  int ierr = 0;         
  TrackStateJacobian Jinv;
  ierr = TrackModel->TrackJacobian( xref, Surf, nSurf, Jinv);  
  if (ierr != 0) {
    return -1;
  }	
       
  // Add scatter noise
  // -----------------------------------
        
  // Local Scatter gain matrix  
  TrackStateGain Gl = TrackModel->GetScatterGain(nxref);   // Backward form -> use nxref not xref
                 
         
  x0 = (Jinv*x0).eval();
  C0.noalias() = (Jinv * C0).eval() *Jinv.transpose() + theta2 * Gl*Gl.transpose();
               
  return 0; 
}
 

/** Returns the predicted chi2 for hit
 */
double TBKalmanB::GetPredictedChi2(const TrackState& p, const TrackStateCovariance& C, const TBHit& hit)
{
  // Get measurement
  const Vector2d& m = hit.GetCoord();
  const Matrix2d& V = hit.GetCov();
  // Compute residual 
  const StateHitProjector& H = GetHMatrix();
  Vector2d r = m - H*p;
  // Compute chisq     	    
  return GetPredictedChi2(r,H,C,V);
}
  
 
/** Returns the chi2 increment for hit
 */
double TBKalmanB::GetChi2Increment(const TrackState& p, const TrackStateCovariance& C, const TBHit& hit)
{
  // Get measurement
  const Vector2d& m = hit.GetCoord();
  const Matrix2d& V = hit.GetCov();
  // Compute residual 
  const StateHitProjector& H = GetHMatrix();
  Vector2d r = m - H*p;
  // Compute chisq     	    
  return GetChi2Increment(r,H,C,V); 
}

/** Returns the predicted chi2
 */
double TBKalmanB::GetPredictedChi2(const Vector2d& r, const StateHitProjector& H, const TrackStateCovariance& C, const Matrix2d& V)
{
  // Residuals weight: W=(V + HCH^T)^-1
  bool invertible = true;
  Matrix2d W=Matrix2d::Zero();// = (V + H*C*H.transpose()).inverse();
  (V + H*C*H.transpose()).computeInverseWithCheck(W,invertible);  // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
  if (!invertible) {
    streamlog_out(MESSAGE3) << "CHi2Increment: matrix inversion failed"
                         << std::endl;
    return -1;
  }	
  // chisq= r^T*W*r
  return (r.transpose()*W*r)[0];
}


/** Returns the chi2 increment
 */
double TBKalmanB::GetChi2Increment(const Vector2d& r, const StateHitProjector& H, const TrackStateCovariance& C, const Matrix2d& V)
{
  // Residuals weight: W=(V - HCH^T)^-1
  bool invertible = true;
  Matrix2d W =Matrix2d::Zero();//= (V - H*C*H.transpose()).inverse();
  (V - H*C*H.transpose()).computeInverseWithCheck(W,invertible);  // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
  if (!invertible) {
    streamlog_out(MESSAGE3) << "CHi2Increment: matrix inversion failed"
                         << std::endl;
    return -1;
  }	
  // chisq= r^T*W*r
  return (r.transpose()*W*r)[0];
}

/** Set number of degrees of freedom in trk
 */
void TBKalmanB::SetNdof(TBTrack& trk)
{

  int ndof = 0;
  int nhits = trk.GetNumHits();

  if ( !TrackModel->BfieldON() ) {  
    if(nhits>2) ndof = 2*nhits-4;
  } else {
    if(nhits>3) ndof = 2*nhits-5;
  }

  trk.SetNDF( ndof );
  
}

/** Compute a priori predicted estimate on first sensor. Returns reference state at first sensor.
 *
 * This function overrides the initial track state parameters and its covariance matrix
 */
TrackState TBKalmanB::ComputeBeamConstraint(ReferenceFrame& FirstSensorFrame, double mom, double charge)
{
   
  // First, we must construct a beam uvw plane in front of the first sensor where the beam 
  // constrained is defined. But this 10 mm in front of the first sensor. 
  ReferenceFrame BeamFrame;
  
  Vector3d  BeamPosition;
  BeamPosition << 0, 0, FirstSensorFrame.GetZPosition()-10;  
  BeamFrame.SetPosition(BeamPosition); 
    
  Matrix3d BeamRotation;
  FillRotMatrixKarimaki(BeamRotation, 0,0,0);
  BeamFrame.SetRotation(BeamRotation); 
  
  // Construct a Gaussian beam state     
  TrackState x_beam;
  x_beam << 0, 0, 0, 0, charge/mom;
  
  // Extrapoate beam  state to first sensor 
  bool error = false;
  TrackState x_first = TrackModel->Extrapolate(x_beam, BeamFrame, FirstSensorFrame, error);  
      
  return x_first;
}


} // Namespace;
