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
  
  ndim = 5; // dimension of state vector
   
  // Project hit coord out of track state
  H = StateHitProjector::Zero();
  H(0,2) = 1.;
  H(1,3) = 1.;

  // Set the initial track state vector
  x0 = TrackState::Zero();
    
  // Set the initial track state covariance matrix
  C0 = TrackStateCovariance::Zero();
    
  for (int i = 0; i < 2; i++) {
    // There is a tradeoff between choosing large values (which makes the
    // tracking numerically unstable) and choosing small values (which biases
    // the track from the seed). I found 1E4 working for this model, but this
    // might need to get tuned.
    C0(i,i) = 1E-2;
    C0(i+2,i+2) = 1E2; 
  } 
  C0(4,4) = 1;  
  
  // Electron mass, in GeV 
  mass = 0.000511;

  // Electron charge in, e
  charge = -1;
  
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
    TrackState Next_State = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);  
      
    if (error) {
      TE.SetCrossed(false);
    } else {
      TE.SetCrossed(true); 
      TE.GetState().GetPars() = Next_State;
    }
    
  } 

  return false;
}


/** Performs track fitting. Returns error flag.
 */
bool TBKalmanB::Fit(TBTrack& trk)
{
   
  // Before the actual fitting, some things need clarification:
  // 
  // 1) The track state estimate (x, C) consists of an estimated parameters  
  //    x and a covariance matrix C.
  // 2) The track state is local, i.e. it is defined with respect to the 
  //    w=0 surface local sensor frame. 
  
  double chisqu = -1;    
  
  // Particle hypothesis  
  mass = trk.GetMass();
  charge = trk.GetCharge();
  
  for(int iter=0; iter<NumIt; iter++) { 
    
    // Get paramters of reference trajectory   
    TrackState RTrkState = trk.GetReferenceState().GetPars();      
    ReferenceFrame RTrkSurf = trk.GetTE( trk.GetReferenceState().GetPlane() ).GetDet().GetNominal();
    
    if ( m_useBC ) {
      // The idea is to use the beam constraint as a reference
      // trajectory 
      double mom = trk.GetMomentum();    
      RTrkSurf = trk.GetTE(0).GetDet().GetNominal();
      RTrkState = ComputeBeamConstraint( x0, C0, RTrkSurf, mass, mom, charge);
    }
    
    // Compute intersections of the reference trajectory
    // with sub detectors. 
     
    // List of crossed track elements  
    vector<int> CrossedTEs;  
    
    // Reference trajectory  
    vector<TrackState> RefStateVec;
    
    // Loop over all track elements 
    int nTE = trk.GetNumTEs();  
    
    for(int iTE=0; iTE<nTE; ++iTE) {
      
      // Next surface along beam line  
      const ReferenceFrame& Next_Surf = trk.GetTE(iTE).GetDet().GetNominal();
           
      // Extrapoate reference state 
      bool error = false;
      auto Next_State = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);  
      
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
    int nCrossed = (int) CrossedTEs.size(); 
    
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
    
    FilterStateVec FFilterData(nCrossed);
    FFilterData[0].Pr_x = x0;
    FFilterData[0].Pr_C = C0;
    
    // Forward FILTER -----------------------------------------------------
    chisqu = FilterPass(trk, CrossedTEs, RefStateVec, +1, FFilterData); 

    // A very basic consistency test 
    if ( std::isnan(chisqu) ||  chisqu < 0 )  {
      trk.SetChiSqu(-1);
      SetNdof(trk);
      return true;
    }

    
    // Init backward filter 
    // The backward pass is initialized using 
    // the smoothed track state at the last 
    // sensor from forward pass. 
    
    FilterStateVec BFilterData(nCrossed);
    BFilterData[nCrossed-1].Pr_x = x0;   
    BFilterData[nCrossed-1].Pr_C = C0;  
    
    // Backward FILTER  -----------------------------------------------------
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
      
      // Get track element 
      TBTrackElement& TE = trk.GetTE( iTE ); 
      
      // Compute unbiased smoothed estimate 
      // of local track state. 
          
      auto& xs =TE.GetState().GetPars(); 
      auto& Cs = TE.GetState().GetCov();
      
      //cout << " plane " << is << " Cf is " << FFilterData[is].Pr_C << endl;
      //xs = FFilterData[is].Pr_x;    
      //Cs = FFilterData[is].Pr_C;

      
      // Smoothing 
      if (is == 0) {
          xs = BFilterData[is].Pr_x; 
          Cs = BFilterData[is].Pr_C;  
      } else if (is == nCrossed-1) {
          xs = FFilterData[is].Pr_x;    
          Cs = FFilterData[is].Pr_C;
      } else {
          auto xb = BFilterData[is].Pr_x; 
          auto Cb = BFilterData[is].Pr_C;
          auto xf = FFilterData[is].Pr_x;    
          auto Cf = FFilterData[is].Pr_C;
                
          bool error = GetSmoothedData(xb, Cb, xf, Cf, xs, Cs);
          if ( error ) {
            trk.SetChiSqu(-1);
            SetNdof(trk); 
            return true; 
          }
                     
      }
      
      // Linearization: This adds the reference parameters         
      xs += RefStateVec[is];
      
      // Compute local ChiSqu   
      if ( TE.HasHit() ) {     
        TBHit& hit =  TE.GetHit();            
        TE.SetChiSqu(GetPredictedChi2(xs, Cs, hit));
      } 
    }  
     
    // Linearization 
    // 
    // We can exploit the current fit to improve the linearization point 
    // iteratively. Note that the beam constraint overwrites the new 
    // linearization point. 
         
    trk.GetReferenceState().Pars = trk.GetTE( CrossedTEs[0] ).GetState().GetPars();  
    trk.GetReferenceState().SetPlane(CrossedTEs[0]);  
    trk.SetMomentum(std::abs(charge/trk.GetTE( CrossedTEs[0] ).GetState().GetPars()[4]));
     
  } // end iterations
     
  // Everything is ok
  trk.SetChiSqu(chisqu);
  SetNdof(trk);  
  return false;
}

/** Filters a new hit 
 */
double TBKalmanB::FilterHit(TBHit& hit, TrackState& xref, TrackState& x0, TrackStateCovariance& C0) 
{ 
  
  // Here, the measurment update takes place  
  double predchi2 = 0; 
    
  // Measured hit coordinates, 2x1 matrix 
  auto& m = hit.GetCoord();
        
  // Covariance for hit coordinates, 2x2 matrix 
  auto& V = hit.GetCov();
  
  bool invertible;
  Matrix2d W ;
  (V + H*C0*H.transpose()).computeInverseWithCheck(W,invertible);  // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
  if (!invertible) {
    streamlog_out(ERROR) << "Hit filtering: Matrix inversion failed. Quit fitting!"
                         << std::endl;
    return -1;
  }	
     
  // This is the predicted residual 
  auto r = m - H*x0 - H*xref; 
      
  // This is the predicted chi2 
  auto chi2mat = r.transpose()*W*r;
  predchi2 = chi2mat[0];   
      
  // Kalman gain matrix K 
  auto K = C0 * H.transpose() * W; 
       
  // This is the filtered state
  x0 += K * r;
  C0 -= C0*H.transpose()*W*H*C0.transpose(); 
        
  return predchi2;
}


/** Propagte state from trackelement te to next track element nte 
 */
int TBKalmanB::PropagateState(TBTrackElement& te, TBTrackElement& nte, TrackState& xref, TrackState& nxref, TrackState& x0, TrackStateCovariance& C0)
{ 
       
  // Reference frame for next track element          
  ReferenceFrame& nSurf = nte.GetDet().GetNominal();

  // Reference frame for current track element          
  ReferenceFrame& Surf = te.GetDet().GetNominal();
  
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
    ierr = MAP_FORWARD( theta2_det, xref, Surf, xref_air, Surf_air, x0, C0 );
    if (ierr != 0) {
      streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
      return -1;
    }	
         
    // MAP estimate [x0,C0] from air to next te 
    // ---------------------------------------
    double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge);   
    ierr = MAP_FORWARD( theta2_air, xref_air, Surf_air, nxref, nSurf, x0, C0 );
    if (ierr != 0) {
      streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
      return -1;
    }	
          
  } else {
        
    // MAP estimate [x0,C0] from old det to air surface
    // ---------------------------------------
    double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge) ;   // Backward form 
    ierr = MAP_BACKWARD( theta2_air, xref, Surf, xref_air, Surf_air, x0, C0 );
    if (ierr != 0) {
      streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
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
      streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
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
    ReferenceFrame& Surf = te.GetDet().GetNominal();      
    
    // Get current reference state xref
    auto& xref = RefStateVec[is];    

    // Copy predicted [x,C]
    auto x0 = Result[is].Pr_x;  
    auto C0 = Result[is].Pr_C;
     
    // Here, the measurment update takes place  
    double predchi2 = 0; 
    
    if ( te.HasHit() ) {
      
      // Measured hit coordinates, 2x1 matrix 
      auto m = te.GetHit().GetCoord();
        
      // Covariance for hit coordinates, 2x2 matrix 
      auto& V = te.GetHit().GetCov();
      
      // Linearization: Measured deviation from reference 
      m -= H*xref; 
            
      // Weigth matrix of measurment 
      bool invertible;
      Matrix2d W ;
      (V + H*C0*H.transpose()).computeInverseWithCheck(W,invertible);  // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
      if (!invertible) {
        streamlog_out(ERROR) << "Hit filtering: Matrix inversion failed. Quit fitting!"
                             << std::endl;
        return -1;
      }	
       
      // This is the predicted residual 
      auto r = m - H*x0; 
      
      // This is the predicted chi2 
      auto chi2mat = r.transpose()*W*r;
      predchi2 = chi2mat[0];   
      
      // Kalman gain matrix K 
      auto K = C0 * H.transpose() * W; 
       
      // This is the filtered state
      x0 += K * r;
      C0 -= C0*H.transpose()*W*H*C0.transpose();        
    }  
    
    // Store results and update chisqu
    chisqu += predchi2;
    
    
     
    // Propagate to next surface
    int inext = is+idir;
  
    if (inext!=istop+idir)  {
         
      // Next surface along filter direction 
      TBTrackElement& nte = TEVec[ CrossedTEs[inext] ]; 
      
      // Parameters for next surface          
      ReferenceFrame& nSurf = nte.GetDet().GetNominal();

      // Get reference state  
      auto& nxref = RefStateVec[inext];  

      // This fitter takes into account scatter in air
      // gaps between sensors. Therefor, we add a virtual 
      // air surface between the detectors which scatters 
      // the track with a material budget equal to the 
      // extrapolation step length between the sensors.
        
      // Extrapolate to air surface between detector planes 
      auto xref_air = xref; 
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
        ierr = MAP_FORWARD( theta2_det, xref, Surf, xref_air, Surf_air, x0, C0 );
        if (ierr != 0) {
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
         
        // MAP estimate [x0,C0] from air to new det surface
        // ---------------------------------------
        double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge);   
        ierr = MAP_FORWARD( theta2_air, xref_air, Surf_air, nxref, nSurf, x0, C0 );
        if (ierr != 0) {
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
        
        
      } else {
        
        // MAP estimate [x0,C0] from old det to air surface
        // ---------------------------------------
        double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge) ;   // Backward form 
        ierr = MAP_BACKWARD( theta2_air, xref, Surf, xref_air, Surf_air, x0, C0 );
        if (ierr != 0) {
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
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
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
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
bool TBKalmanB::GetSmoothedData( TrackState& xb, TrackStateCovariance& Cb, TrackState& xf, TrackStateCovariance& Cf,
                                TrackState& Smoothed_State, TrackStateCovariance& Smoothed_Cov)
{ 
  
  // Error flag  
  bool invertible;
  
  // Compute forward weight matrix
  TrackStateWeight fW; 
  Cf.computeInverseWithCheck(fW,invertible);  
  if (!invertible) {
    streamlog_out(ERROR) << "Smoothing 1: Matrix inversion failed. Quit fitting!"
                         << std::endl;
    return invertible;
  }	
  
  // Compute backward weight matrix
  TrackStateWeight bW;
  Cb.computeInverseWithCheck(bW,invertible);  
  if (!invertible) {
    streamlog_out(ERROR) << "Smoothing 2: Matrix inversion failed. Quit fitting!"
                         << std::endl;
    return invertible; 
  }	 
  
  // Weighted (smoothed) covariance   
  (fW + bW).computeInverseWithCheck(Smoothed_Cov,invertible);  
  if (!invertible) {
    streamlog_out(ERROR) << "Smoothing 3: Matrix inversion failed. Quit fitting!"
                         << std::endl;
    return invertible;
  }	
  
  // Weighted state           
  Smoothed_State = Smoothed_Cov*(fW*xf + bW*xb);
  
  // Successfull return  
  return invertible;  
}




int TBKalmanB::MAP_FORWARD(  double theta2, 
                        TrackState& xref, ReferenceFrame& Surf,
                        TrackState& nxref, ReferenceFrame& nSurf,
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
  TrackStateJacobian J;      
  TrackModel->TrackJacobian( xref, Surf, nSurf, J);  
  
  auto C1 = J*C0*J.transpose(); 
        
  // Add scatter noise
  // -----------------------------------
        
  // Local Scatter gain matrix      
  TrackStateGain Gl = TrackModel->GetScatterGain(xref);

  //cout << "BENNI HACK scatter gain matrix " << Gl << endl; 
            
  // General Scatter gain matrix 
  auto G = J*Gl;   
  
  // Variance of projected scatter angles   
  auto Q=theta2*Matrix2d::Identity();
             
  C1 += G*Q*G.transpose();     
         
  x0 = J*x0;  
  C0 = C1;  

  
          
  return 0; 
}

int TBKalmanB::MAP_BACKWARD(  double theta2, 
                        TrackState& xref, ReferenceFrame& Surf,
                        TrackState& nxref, ReferenceFrame& nSurf,
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
  TrackStateJacobian Jinv;
  TrackModel->TrackJacobian( xref, Surf, nSurf, Jinv);  
  
  auto C1 = Jinv * C0 *Jinv.transpose();
        
  // Add scatter noise
  // -----------------------------------
        
  // Local Scatter gain matrix  
  TrackStateGain Gl = TrackModel->GetScatterGain(nxref);   // Backward form -> use nxref not xref
                 
  // Variance of projected scatter angles   
  auto Q=theta2*Matrix2d::Identity();

             
  C1 +=  Gl*Q*Gl.transpose();    // Backward form -> use Gl not G
         
  x0 = Jinv*x0;  
  C0 = C1;  
               
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
double TBKalmanB::GetChi2Increment(TrackState& p, TrackStateCovariance& C, TBHit& hit)
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
  bool invertible;
  Matrix2d W ;
  (V + H*C*H.transpose()).computeInverseWithCheck(W,invertible);  // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
  if (!invertible) {
    streamlog_out(ERROR) << "CHi2Increment: matrix inversion failed"
                         << std::endl;
    return -1;
  }	
  // chisq= r^T*W*r
  auto chisq = r.transpose()*W*r;
  return chisq[0];
}


/** Returns the chi2 increment
 */
double TBKalmanB::GetChi2Increment(const Vector2d& r, const StateHitProjector& H, const TrackStateCovariance& C, const Matrix2d& V)
{
  // Residuals weight: W=(V - HCH^T)^-1
  bool invertible;
  Matrix2d W ;
  (V - H*C*H.transpose()).computeInverseWithCheck(W,invertible);  // HCH^T is only one 2x2 block from C if H is a simple projectior. That could be done better i guess.
  if (!invertible) {
    streamlog_out(ERROR) << "CHi2Increment: matrix inversion failed"
                         << std::endl;
    return -1;
  }	
  // chisq= r^T*W*r
  auto chisq = r.transpose()*W*r;
  return chisq[0];
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
 */
TrackState TBKalmanB::ComputeBeamConstraint( TrackState& x0, TrackStateCovariance& C0, ReferenceFrame& FirstSensorFrame, double mass, double mom, double charge)
{
   
  // First, we must construct a beam uvw plane in front of the first sensor where the beam 
  // constrained is defined. 
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
  auto x_first = TrackModel->Extrapolate(x_beam, BeamFrame, FirstSensorFrame, error);  
  
  // MAP estimate [x_beam,C_beam] from beam frame to first sensor
  int ierr = 0; 
  double length = TrackModel->GetSignedStepLength(x_beam, BeamFrame, FirstSensorFrame); 
  double theta2_air = materialeffect::GetScatterTheta2(x_beam, length, materialeffect::X0_air, mass, charge);   
  ierr = MAP_FORWARD( theta2_air, x_beam, BeamFrame, x_first, FirstSensorFrame, x0, C0 );
  if (ierr != 0) {
    streamlog_out(ERROR) << "ERR: Problem with track extrapolation"
                         << std::endl;
  }	
      
  return x_first;
}


} // Namespace;

