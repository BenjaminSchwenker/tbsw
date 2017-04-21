// TBKalmanSRIF implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Marlin includes 
#include <marlin/Global.h>
#include <streamlog/streamlog.h>
  
// DEPFETTrackTool includes
#include "TBKalmanSRIF.h"
#include "MaterialEffect.h"
#include "StraightLineTrackModel.h"

#include "ThreeDModel.h" 

// C++ standard library includes 
#include <cmath>
#include <limits>

// Namespaces
using namespace CLHEP;
using namespace marlin;
using namespace std; 

namespace depfet {

/** Constructor 
 */
TBKalmanSRIF::TBKalmanSRIF()
{ 
  // Number of iterations for Kalman filter
  NumIt = 1;     
  
  // By default, no cut is used
  OutlierCut = numeric_limits< double >::max();          
  
  // By default, beam constraint is not used 
  m_useBC = false; 
  
  m_sizex = 0;       // spot size, mm
  m_sizey = 0;       // spot size, mm
  m_divx = 0;     // divergence, rad
  m_divy = 0;     // divergence, rad
  m_corrx = 0;      // beam correlation coeff  
  m_corry = 0;      // beam correlation coeff
  
  // Project hit coord out of track state
  H = HepMatrix(2,4,0);
  H[0][2] = 1.;
  H[1][3] = 1.; 
  
  TrackModel = new StraightLineTrackModel();   
}

TBKalmanSRIF::~TBKalmanSRIF()
{ 
  delete TrackModel; 
}

/** Extrapolates the track seed to all planes. Returns error flag. 
*/
bool TBKalmanSRIF::ExtrapolateSeed(TBTrack& trk) 
{
 
  // Get seed parameters   
  HepMatrix RTrkState = trk.GetReferenceState().GetPars();      
  ReferenceFrame RTrkSurf;
    
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
    HepMatrix Next_State = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);
    if (error) {
        return error;   
    }
             
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
bool TBKalmanSRIF::Fit(TBTrack& trk)
{
   
  // Before the actual fitting, some things need clarification:
  // 
  // 1) The Kalman filter is linearized around a reference trajectory.
  //    So, it 'filters' not the track state x, but the deviation  
  //    from a reference track state.  
  // 1a)We will call the deviation from the reference trajectory also 
  //    a (track) state. Do not forget to add reference state later.  
  // 1b)The track state estimate (x, C) consists of an estimated parameters  
  //    x and a covariance matrix C. 
  // 1c)All scatter gain and state transition matrices are calculated for
  //    the reference track states.
  // 1d)Fit iterations means to carry out a fit with a better reference 
  //    trajectory, or after removing outlier measurments.  
  // 2) The Kalman filter is mechanized as a square root information
  //    filter (SRIF). A SRIF is not directly using variables [x,C], 
  //    but works instead with 'data equations' of form z=R*x+v. 
  // 2a)Here: R is a square root of the information matrix Y=1/C. 
  //    That is Y=R^tR. We will call R data matrix.
  // 2b)The data vector z is defined as z=Rx. It plays a similiar role 
  //    is an observation vector in classical least square fitting.
  // 2c)The error v is a white noise with unit variance and zero mean. 
  //    Note that z is dimensionless.   
  // 2d)Going from (x,C) to (z,R) normalizes all state variables.
  //    This tends to reduce the dynamic range of the variables in all 
  //    computations. 
  // 
  // Reference: 
  // The SRIF method is very nicely explained in the Phd thesis of Paul 
  // Kaminsky. Compared to a direct implementation of Kalman's matrix
  // equation ('conventional filter'), the SRIF has order of magnitude
  // better numerical stability and accuracy. 
   
  double chisqu = -1;    
  trk.SetChiSqu(-1);
  
  for(int iter=0; iter<NumIt; iter++) { 
    
    //cout << "iter " << iter << endl; 
    
    // There are only two things that change between interations:
    //  
    // a) Hits in the track using TBTrackElement.RemoveHit()
    // b) Track reference state using SetReferenceState() 

    // Initialize Kalman filters 
    // 
    // We need to initialiaze the forward and backward filters.
    // At the beginning, we use z=0 and R=0. This is equivalent 
    // to using a very large initial covariance C0 -> infinity. 
    // This eliminates any bias of finite states to the initial 
    // track parameter estimate 
    
    HepMatrix z0f(4,1,0);
    HepMatrix R0f(4,4,0);
    
    HepMatrix z0b(4,1,0);
    HepMatrix R0b(4,4,0);
      
    // Idea is to use reference trajectory in TBTrack object from the 
    // pattern reco as linearization point in fitter
      
    HepMatrix RTrkState = trk.GetReferenceState().GetPars();         
    ReferenceFrame RTrkSurf; 
    
    if ( m_useBC ) {
      // Idea is to use the beam constraint as a reference
      // trajectory 
      double mass = trk.GetMass();
      double mom = trk.GetMomentum(); 
      int charge = trk.GetCharge();
      
      RTrkSurf = trk.GetTE(0).GetDet().GetNominal();
      RTrkState = ComputeBeamConstraint( z0f, R0f, RTrkSurf, mass, mom, charge);
    }
    
    // Compute intersections of the reference trajectory
    // with sub detectors. 
    
    // List of crossed track elements  
    vector<int> CrossedTEs;  
    
    // Reference trajectory  
    vector<HepMatrix> RefStateVec;
    
    // Loop over all track elements 
    int nTE = trk.GetNumTEs(); 
      
    for(int iTE=0; iTE<nTE; ++iTE) {
       
      // Next surface along beam line  
      const ReferenceFrame& Next_Surf = trk.GetTE(iTE).GetDet().GetNominal();
      
      // Extrapoate reference state 
      bool error = false;
      HepMatrix Next_State = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);
      if (error) {
        return error;   
      }
      
      // The detector plane is crossed by the particle iff there is either a hit 
      // on the detector or the reference track crosses the module area.  
      double u = Next_State[2][0]; 
      double v = Next_State[3][0];
      
      if ( trk.GetTE(iTE).HasHit() || trk.GetTE(iTE).GetDet().ModuleCrossed(u,v) ) {
        trk.GetTE(iTE).SetCrossed(true);
        // Remember plane number 
        CrossedTEs.push_back(iTE);
        // Remember reference state 
        RefStateVec.push_back(Next_State);      
      } else {
        trk.GetTE(iTE).SetCrossed(false);
      }
          
    } 
     
    // Number of crossed track elements  
    int nCrossed = (int) CrossedTEs.size(); 
   
    // The reference trajectory should at least 
    // intersect with one detector.  
    
    if ( nCrossed == 0 ) {
      trk.SetChiSqu(-1);  
      return true;
    }   
    
    // Init backward filter 
    KalFilterVec BFilterData(nCrossed);
    BFilterData[nCrossed-1].Pr_z = z0b; 
    BFilterData[nCrossed-1].Pr_R = R0b;  

    // Init forward filter 
    KalFilterVec FFilterData(nCrossed);
    FFilterData[0].Pr_z = z0f;
    FFilterData[0].Pr_R = R0f;
     
    // Backward KALMAN FILTER  -----------------------------------------------------
    chisqu = KalmanFilterPass(trk, CrossedTEs, RefStateVec, -1, BFilterData);  
    
    // A very basic consistency test 
    if ( std::isnan(chisqu) ||  chisqu < 0 )  {
      trk.SetChiSqu(-1);
      return true;
    }
    
    // Forward KALMAN FILTER -----------------------------------------------------
    KalmanFilterPass(trk, CrossedTEs, RefStateVec, +1, FFilterData); 
    
    // Kalman Smoothing 
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
         
      HepMatrix zb = BFilterData[is].Pr_z; 
      HepMatrix Rb = BFilterData[is].Pr_R;
      HepMatrix zf = FFilterData[is].Pr_z;    
      HepMatrix Rf = FFilterData[is].Pr_R;
      
      // Get pointers to paramter vector and 
      // covariance matrix
       
      HepMatrix& p0 =TE.GetState().GetPars(); 
      HepSymMatrix& C0 = TE.GetState().GetCov();
              
      // Performs the smoothing
      bool error = GetSmoothedData(zb, Rb, zf, Rf, p0, C0);
      
      // Adds reference parameters         
      p0 += RefStateVec[is];
      
      if ( error ) {
        streamlog_out(MESSAGE2) << "TRK ERROR: Smoothing has failed at track element " << iTE << std::endl;
        TE.SetChiSqu(0);
      } 
              
      // Compute local ChiSqu   
      if ( TE.HasHit() && !error ) {     
        TBHit& hit =  TE.GetHit();            
        TE.SetChiSqu(GetPredictedChi2(p0, C0, hit));
      } 
         
    } 
           
    // Outlier logic
    // 
    // Remove hit with largest smoother chi2 value, i.e.
    // the worst measurment. This requires one more 
    // fitter iteration!
    
    double WorstChi2 = 0;   // Worst Chisqu  
    int WorstTE = -1;       // Worst TE
    
    for(int is=0; is<nCrossed; ++is) {
       
      // Get plane number 
      int iTE = CrossedTEs[is];  
      
      // Get track element 
      TBTrackElement& TE = trk.GetTE( iTE ); 
     
      // Remember TE with max chi2
      if ( WorstChi2 < TE.GetChiSqu() ) {
        WorstChi2 = TE.GetChiSqu(); 
        WorstTE = iTE;   
      }
    }
    
    // One iteration is needed to refit
    if ( iter < NumIt-1 &&  WorstChi2 > OutlierCut ) { 
      // Exclude worst hit from track   
      trk.GetTE( WorstTE ).RemoveHit();
    }
    
    // Linearization point
    // 
    // We can exploit the current fit to improve the linearization point 
    // iteratively. Note that the beam constraint overwrites the new 
    // linearization point. 
         
    HepMatrix TmpState = trk.GetTE( CrossedTEs[0] ).GetState().GetPars();       
    ReferenceFrame TmpSurf = trk.GetTE( CrossedTEs[0] ).GetDet().GetNominal(); 
    
    // Extrapolate State to Z=0 
    bool error;
    HepMatrix NewRef = TrackModel->Extrapolate(TmpState, TmpSurf, ReferenceFrame(), error);
    if (error) {
        streamlog_out(MESSAGE2) << "ERROR: Propagation to next detector failed. Quit fitting!"
                                << std::endl;
        // Quit fitter, set error flag
        trk.SetChiSqu(-1);
        return error;       
    }
    
    trk.GetReferenceState().GetPars() = NewRef;     
     
  } // End fit iterations 

  // Everything is ok
  trk.SetChiSqu(chisqu);
  return false;
}


/** Run Kalman filter on track. Returns fit chi2. 
 *    
 *  Input: StateVec    : Vector of reference track parameters  
 *  Input: IDIR        : -1 for backward filter, 1 for forward filter 
 *  Output : RESULT    : Estimated [R,z] pairs for all detectors 
 */
double TBKalmanSRIF::KalmanFilterPass(TBTrack& trk, std::vector<int>& CrossedTEs, REFTrack& StateVec, int idir, KalFilterVec& Result) 
{ 
  
  // To start ierr is set to 0 (= OK)
  int ierr = 0; 
  // Final chisqu of track fit
  double chisqu = 0;
       
  // Particle hypothesis 
  double mom = trk.GetMomentum(); 
  double mass = trk.GetMass();
  int charge = trk.GetCharge();
  
  // Process vector of track elements
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
    HepMatrix& xref = StateVec[is];    
    
    // Copy predicted [R0,z0]
    HepMatrix z0 = Result[is].Pr_z;  
    HepMatrix R0 = Result[is].Pr_R;
     
    //cout << "Predicted z0 " << z0 << endl;
    //cout << "Predicted R0 " << R0 << endl;
      
    // MEASUREMENT UPDATE 
     
    double predchi2 = 0; 
    
    if ( te.HasHit() ) {  

      // Measured hit coordinates, 2x1 matrix 
      HepMatrix m = te.GetHit().GetCoord();
      
      // Covariance for hit coordinates, 2x2 matrix 
      HepSymMatrix& V = te.GetHit().GetCov();
       
      // Measured deviation from reference 
      // trajectory 
      HepMatrix zm = m - H*xref; 
        
      // Measurement tableau TAB represents an 
      // augmented data equation including the 
      // measured hits. 
      // 
      //  TAB:
      //  
      //  ( R0 | z0 )  
      //  ( H  | zm )     
      // 
      // Note that the data equation is not unique. 
      // Data equations can be transformed into each
      // other using orthogonal transformations. We
      // can exploit this freedom to solve the data 
      // equation partially for zm 
      // 
      //  ( R0' | z0' )  
      //  ( 0   | e   ) 
      // 
      // In particular, this gives a new data equation
      // for track state plus a measurement error e.   
      
      HepMatrix TAB(6,5,0);
      
      // Copy R0 and z0 into TAB
      for (int l=0;l<4;++l) {
        TAB[l][4] = z0[l][0]; 
        for (int m=0;m<4;++m) {
          TAB[l][m] = R0[l][m]; 
        }
      }
      
      // Add u measurment - rescale to unit variance 
      TAB[4][4] = zm[0][0]/std::sqrt(V[0][0]); 
      TAB[4][2] = 1./std::sqrt(V[0][0]); 
       
      // Add v measurment - rescale to unit variance 
      TAB[5][4] = zm[1][0]/std::sqrt(V[1][1]); 
      TAB[5][3] = 1./std::sqrt(V[1][1]);
     
      // Carry out Householder tridiagonalization on TAB.
      //
      // From CLHEP documentation: 
      // void house_with_update(HepMatrix*A,int row=1, int col=1)
      // 
      // It changes A to PA, where P is the Householder trafo to 
      // zero all elements of A in column col from row+1 last row
      // of A. 
      // NOTE: Columns/rows are indexed starting at 1!!!!
      // NOTE: w is a column vector defining P = 1-beta*w*wT
      // ERROR: CLHEP Householder crashes, if w==0. Needs startup!
       
      // The idea is to zero the 4 columns below the diagonal of TAB
      // using Housholder transformations. 
      
      //cout << "initial TAB " << TAB << endl;
      
      int startcol = StartupHouseholder(TAB, 4); 
      for (int icol=startcol;icol<=4;++icol){ 
        house_with_update(&TAB,icol,icol); 
        //cout << "TAB " << TAB << endl;
      }
      
      //cout << "final TAB " << TAB << endl; 
       
      // Actually overwrite data equation
      for (int l=0;l<4;++l) {
        z0[l][0] = TAB[l][4];  
        for (int m=0;m<4;++m) {
          R0[l][m] = TAB[l][m]; 
        }
      }
       
      // And the predicted chi2 is for free :) 
      predchi2 = TAB[4][4]*TAB[4][4] + TAB[5][4]*TAB[5][4]; 
    }  
    
    //cout << "Filtered chisq " << predchi2 << endl;
    //cout << "Filter z0 " << z0 << endl;
    //cout << "Filter R0 " << R0 << endl;
     
    // Store results and update chisqu
    chisqu += predchi2;
    Result[is].Chi2Pred = predchi2; 
    Result[is].Up_z = z0;  
    Result[is].Up_R = R0;
    
    // TIME UPDATE  
    // 
    // Estimated track parameters are mapped to the 
    // next crossed detector in filter direction. 
    // 
    // The fitter takes into account scatter in air  
    // between sensors. Therefor, we add a virtual 
    // air surface between the detectors which scatters 
    // the track with a material budget equal to the 
    // extrapolation step length between the sensors.   
     
    // Index of next crossed detector   
    int inext = is+idir;
      
    if (inext!=istop+idir)  {
        
      // Next crossed TE along filter direction 
      TBTrackElement& nte = TEVec[ CrossedTEs[inext] ]; 
      
      // Get surface parameters          
      ReferenceFrame& nSurf = nte.GetDet().GetNominal();
      
      // Get reference state  
      HepMatrix& nxref = StateVec[inext];    
      
      // Get signed flight length in air between detectors 
      double length = TrackModel->GetSignedStepLength(xref, Surf, nSurf);
      
      // Compute air surface and reference state 
      HepMatrix xref_air = xref; 
      ReferenceFrame Surf_air = Surf; 
      TrackModel->Extrapolate(xref_air, Surf_air, length/2);
      
      if ( idir > 0 ) {  
      
        // MAP estimate [z0,R0] from old det to air surface
        // ---------------------------------------
        double l0 = te.GetDet().GetTrackLength(xref[2][0], xref[3][0], xref[0][0], xref[1][0]);
        double X0 = te.GetDet().GetRadLength(xref[2][0],xref[3][0]);    
        double theta2_det = materialeffect::GetScatterTheta2(xref, l0, X0, mass, charge);   
        ierr = MAP_FORWARD( theta2_det, xref, Surf, xref_air, Surf_air, z0, R0 );
        if (ierr != 0) {
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
         
        // MAP estimate [z0,R0] from air to new det surface
        // ---------------------------------------
        double theta2_air = materialeffect::GetScatterTheta2(xref_air, length, materialeffect::X0_air, mass, charge);   
        ierr = MAP_FORWARD( theta2_air, xref_air, Surf_air, nxref, nSurf, z0, R0 );
        if (ierr != 0) {
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
        
      } else {
        
        // MAP estimate [z0,R0] from old det to air surface
        // ---------------------------------------
        double theta2_air = materialeffect::GetScatterTheta2(xref, length, materialeffect::X0_air, mass, charge) ;   // Backward form 
        ierr = MAP_BACKWARD( theta2_air, xref, Surf, xref_air, Surf_air, z0, R0 );
        if (ierr != 0) {
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
        
        // MAP estimate [z0,R0] from air to new det surface
        // ---------------------------------------
        double l0 = nte.GetDet().GetTrackLength(nxref[2][0], nxref[3][0], nxref[0][0], nxref[1][0]);
        double X0 = nte.GetDet().GetRadLength(nxref[2][0],nxref[3][0]);    
        double theta2_det = materialeffect::GetScatterTheta2( xref_air, l0, X0, mass, charge);   // Backward form    
        ierr = MAP_BACKWARD( theta2_det, xref_air, Surf_air, nxref, nSurf, z0, R0 );          
        if (ierr != 0) {
          streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                               << std::endl;
          return -1;
        }	
        
      }
      
      
      // Store results        
      Result[inext].Chi2Pred = 0; 
      Result[inext].Pr_z = z0;  
      Result[inext].Pr_R = R0;
      
    }        
  } 
    
  // Successfull return  
  return chisqu;  
}

int TBKalmanSRIF::MAP_FORWARD(  double theta2, 
                        HepMatrix& xref, ReferenceFrame& Surf, 
                        HepMatrix& nxref, ReferenceFrame& nSurf, 
                        HepMatrix& z0, HepMatrix& R0
                     )
{
  
  // The linearized mapping from state x_k at k to state x_k+1 at k+1:
  // 
  // x_k+1 = J_k * x_k + G_k * w_k    (*)                    
  // 
  //  J_k     : Transport matrix from k to k+1; lin. around ref state  
  //  w_k     : Vector of scatter angles for scatterers between k->k+1
  //  Q_k     : Diag. covariance matrix for w_k (w_k are zero mean)
  //  G_k     : Scatter gain matrix, wrt. ref state   
 
  int ierr; 
    
  HepMatrix J(4, 4);
  ierr = TrackModel->TrackJacobian( xref, Surf, nSurf, J);
  if (ierr != 0) {
    streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                         << std::endl;
    return -1;
  }	
        
  HepMatrix Jinv(4, 4);
  ierr = TrackModel->TrackJacobian( nxref, nSurf, Surf, Jinv);
  if (ierr != 0) {
    streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                         << std::endl;
    return -1;
  }	 
  
  // Scatter gain matrix  
  HepMatrix Gl(4, 2);    
  TrackModel->GetScatterGain(xref, Gl); 
   
  // For the SRIF implementation, the augmented data equation reads 
  //  
  // MAP:
  //  
  // ( Rw_k     | 0   | 0  )      or      0 = Rw_k*w_k + v_w   
  // ( -Rd*G_k  | Rd  | z_k)            z_k = -Rd*G_k*w_k + Rd*x_k+1 + v_k 
  //  
  // Rd      = R_k * Jinv_k
  // 1/Q_k   = Rw_k^t * Rw_k
  // G_k     = J_k * Gl_k
  //       
  // The first quation represents the a priori about scatterings. The 2nd
  // data equation represents the a priori about track state [R_k,z_k].   
  
  double rw = 1./std::sqrt(theta2);                       
  HepMatrix Rd = R0 * Jinv;     
  HepMatrix G = J*Gl;   
  HepMatrix B = -Rd*G;  
  
  // The augmented data equation fits in tableau MAP 
  HepMatrix MAP(6,7,0);
                                   
  // Fill MAP 
  MAP[0][0] = rw; 
  MAP[1][1] = rw; 
  for (int l=0;l<4;++l) {
    // Fill -G*Rd
    MAP[l+2][0] = B[l][0];
    MAP[l+2][1] = B[l][1];
    // Fill z0
    MAP[l+2][6] = z0[l][0];
    // Fill Rd
    for (int m=0;m<4;++m) {
      MAP[l+2][m+2] = Rd[l][m]; 
    }
  }      
      
  // Again, we will make use of Householder trafos to 
  // put MAP into the more favorable form 
  //  
  // ( Rw(1) | Rwx(1) | zw(1) )  
  // ( 0     | R1     | z1    ) 
      
  int startcol = StartupHouseholder(MAP, 2);
      
  for (int icol=startcol;icol<=2;++icol){ 
    house_with_update(&MAP,icol,icol); 
  }
          
  // Actually overwrite data equation
  for (int l=0;l<4;++l) { 
    z0[l][0] = MAP[l+2][6];
    for (int m=0;m<4;++m) {
      R0[l][m] = MAP[l+2][m+2]; 
    }
  }
     
  return 0; 
}

int TBKalmanSRIF::MAP_BACKWARD(  double theta2, 
                        HepMatrix& xref, ReferenceFrame& Surf, 
                        HepMatrix& nxref, ReferenceFrame& nSurf, 
                        HepMatrix& z0, HepMatrix& R0
                     )
{
  
  // The linearized mapping from state x_k+1 at k+1 to state x_k at k:
  // 
  // x_k = Jinv_k * x_k+1 - Jinv_k * G_k * w_k              (**) 
  // 
  // J_k     : Transport matrix from k to k+1 -> around ref state  
  // w_k     : Vector of all scatter angles for scatterers between k->k+1
  // Q_k     : Diag. covariance matrix for w_k (w_k are zero mean)
  // G_k     : Scatter gain matrix -> around ref state   
  // 
  // Note that the equations are not invariant under time reversal. The scatter
  // angles w_k are defined for a track moving forward in time.  
  
  int ierr; 
  
  HepMatrix J(4, 4);
  ierr = TrackModel->TrackJacobian( nxref, nSurf, Surf, J);
  if (ierr != 0) {
    streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
    return -1;
  }	
        
  HepMatrix Jinv(4, 4);
  ierr = TrackModel->TrackJacobian( xref, Surf, nSurf, Jinv);
  if (ierr != 0) {
    streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                           << std::endl;
    return -1;
  }	 
     
  HepMatrix Gl(4, 2);    
  TrackModel->GetScatterGain(nxref, Gl); // backward SRIF!!! 
    
  // For the SRIF implementation, the augmented data equation reads now  
  //  
  // MAP:
  // 
  // ( Rw_k      | 0       | 0     )      or         0 = Rw_k*w_k + v_w   
  // ( R_k+1*G_k | R_k+1*J | z_k+1 )             z_k+1 = (R_k+1*G_k)*w_k + (R_k+1*J)*x_k + v_k 
  //  
  // 1/Q_k   = Rw_k^t * Rw_k
  // G_k     = J_k * Gl_k
  // 
  // The first quation represents the a priori about scatterings. The 2nd
  // data equation represents the a priori about track state [R_k+1,z_k+1]. 
  
  double rw = 1./std::sqrt(theta2); 
  HepMatrix G = J*Gl;             
  HepMatrix B = R0*G;            // backward SRIF!!!
  HepMatrix A = R0 * J;          // backward SRIF!!!
  
  // The augmented data equation fits in tableau MAP 
  HepMatrix MAP(6,7,0);
     
  // Fill MAP row wise  
  MAP[0][0] = rw; 
  MAP[1][1] = rw; 
  for (int l=0;l<4;++l) {
    // Fill R0*G
    MAP[l+2][0] = B[l][0];
    MAP[l+2][1] = B[l][1];
    // Fill z0
    MAP[l+2][6] = z0[l][0];
    // Fill  R0*J
    for (int m=0;m<4;++m) {
      MAP[l+2][m+2] = A[l][m]; 
    }
  }      
         
  // Again, we will make use of Householder trafos to 
  // put MAP into the more favorable form 
  //  
  // ( Rw(1) | Rwx(1) | zw(1) )  
  // ( 0     | R1     | z1    ) 
      
  int startcol = StartupHouseholder(MAP, 2);
      
  for (int icol=startcol;icol<=2;++icol){ 
    house_with_update(&MAP,icol,icol); 
  }  
      
  // Actually overwrite data equation
  for (int l=0;l<4;++l) { 
    z0[l][0] = MAP[l+2][6];
    for (int m=0;m<4;++m) {
      R0[l][m] = MAP[l+2][m+2]; 
    }
  }
      
  return 0; 
}


/** Returns the predicted chi2 for hit
 */
double TBKalmanSRIF::GetPredictedChi2(CLHEP::HepMatrix& p, CLHEP::HepSymMatrix& C, TBHit& hit)
{
  // Get measurement
  HepMatrix m = hit.GetCoord();  
  HepSymMatrix V = hit.GetCov();
  // Compute residual 
  HepMatrix& H = GetHMatrix();   
  HepMatrix r = m - H*p; 
  // Compute chisq     	    
  return GetPredictedChi2(r,H,C,V); 
}
  
 
/** Returns the chi2 increment for hit
 */
double TBKalmanSRIF::GetChi2Increment(HepMatrix& p, HepSymMatrix& C, TBHit& hit)
{
  // Get measurement
  HepMatrix m = hit.GetCoord();  
  HepSymMatrix V = hit.GetCov();
  // Compute residual 
  HepMatrix& H = GetHMatrix();   
  HepMatrix r = m - H*p; 
  // Compute chisq     	    
  return GetChi2Increment(r,H,C,V); 
}

/** Returns the predicted chi2
 */
double TBKalmanSRIF::GetPredictedChi2( CLHEP::HepMatrix& r, CLHEP::HepMatrix& H, CLHEP::HepSymMatrix& C, CLHEP::HepSymMatrix& V)
{
  // Residuals weight: W=(V - HCH^T)^-1
  int ierr;
  HepSymMatrix W = (V + C.similarity(H)).inverse(ierr);
  if (ierr) {
    streamlog_out(ERROR) << "CHi2Increment: matrix inversion failed"
                         << std::endl;


    return -1;
  }	
  // chisq= r^T*W*r
  HepMatrix chisq = r.T()*W*r; 
  return chisq[0][0];  
}


/** Returns the chi2 increment
 */
double TBKalmanSRIF::GetChi2Increment( HepMatrix& r,
                                HepMatrix& H,
                                HepSymMatrix& C,
                                HepSymMatrix& V)
{ 
  // Residuals weight: W=(V - HCH^T)^-1
  int ierr;
  HepSymMatrix W = (V - C.similarity(H)).inverse(ierr);
  if (ierr) {
    streamlog_out(ERROR) << "CHi2Increment: matrix inversion failed"
                         << std::endl;
    return -1;
  }	
  // chisq= r^T*W*r
  HepMatrix chisq = r.T()*W*r; 
  return chisq[0][0];
}


/** Compute the weighted means of forward and backward filter estimated. In a 
 *  numerically robust square root way. Returns error flag.
 */
bool TBKalmanSRIF::GetSmoothedData( HepMatrix& zb, HepMatrix& Rb, HepMatrix& zf, HepMatrix& Rf,
                                HepMatrix& Smoothed_State, HepSymMatrix& Smoothed_Cov)
{ 
    
  // This is a square root(!) procedure to calculate the 
  // weighted means of the forward filter/backward filter
  // estimates. The problem consists in finding the least 
  // square solution of the following equation  
  // 
  //  WGM:
  //  
  //  ( R_F | z_F )   or  z_F   =  R_F  * x + v_F
  //  ( R_B | z_B )       z_B   =  R_B  * x + v_B
  // 
  // Here, x denotes the true state at current surface, [R_F,z_F] and 
  // [R_B,z_B] are the forward/backward estimates excluding the 
  // hit at the current detector. Then, v_B and v_F denote zero 
  // mean and unit covariance errors of the forward/backward filter. 
  // The task is to compute a least square estimate + covariance for 
  // x given data equations above. 
  // 
  // The square root solution is found using the standard Householder 
  // triangularization procedure. 
  // 
  //  ( R_F  | z_F  )   --> ( R* | z* ) 
  //  ( R_B  | z_B  )       ( 0  | e* )
  // 
  // The final solution is [R*,z*]. Note that R* is always upper triangular. 
  // The the optimal x is just given by back substution from z* = R* x. The 
  // covariance is given by C* = R*^-1 * R*^-t. 
  // 
  // Reference: The idea of a Sqauare Root Two Filter Smoother goes back to 
  // the Phd thesis of Paul Kaminski. Just google it :) 
  
  // Error flag 
  int ierr = 0; 
  
  HepMatrix WGM(8,5,0);
      
  // Copy R's and z's into WGM 
  for (int l=0;l<4;++l) {
    WGM[l][4]   = zf[l][0]; 
    WGM[l+4][4] = zb[l][0];
    
    for (int m=0;m<4;++m) {
      WGM[l][m]   = Rf[l][m]; 
      WGM[l+4][m] = Rb[l][m]; 
    }
  }
      
  //cout << "initial WGM " << WGM << endl;
      
  int startcol = StartupHouseholder(WGM, 4); 
  for (int icol=startcol;icol<=4;++icol){ 
    house_with_update(&WGM,icol,icol); 
    //cout << "WGM " << WGM << endl;
  }
      
  //cout << "final WGM " << WGM << endl; 
       
  // Actually overwrite forward estimate   
  for (int l=0;l<4;++l) {
    zf[l][0] = WGM[l][4];  
    for (int m=0;m<4;++m) {
      Rf[l][m] = WGM[l][m]; 
    }
  }
   
  //cout << "R " << Rf << endl; 
  //cout << "z " << zf << endl; 
  //cout << "determinant(R) " << Rf.determinant() << endl;  
  
  // We solve z = Rx with R upper traingular -> this overwrites z with x 
  back_solve(Rf,&zf);
  //cout << "x " << zf << endl;
   
  // This overwrites R with R^-1 
  Rf.invert(ierr);
  //cout << "Rinv " << R << endl;
  // This overwrites R^-1 with C    
  HepMatrix B =Rf*Rf.T(); 
  //cout << "C " << B << endl;  
  
  // Squaring up A to A*At gives always a symmetric matrix, but CLHEP does not 
  // know this; fucking shit 
  HepSymMatrix ttemp(4,0);
  
  for (int l=0; l<4; ++l){
    for (int m=0; m<=l; ++m) {
      ttemp[l][m] = B[l][m];  
    }
  }  
   
  Smoothed_Cov = ttemp;  
  Smoothed_State = zf; 
   
  // Successfull return  
  return ierr;
     
}

/** Returns a proper start column for Householder triangularization 
 *  
 * The method tries to triangularize first maxcol columns of A 
 * by row permutations. It returns the first column, where this 
 * is not enough.   
 */
int TBKalmanSRIF::StartupHouseholder(CLHEP::HepMatrix& A, int maxcol)
{
  int nrows = A.num_row();
  int ncols = A.num_col();
  
  // loop over columns 
  for (int l=1;l<=maxcol;++l){
    
    int non_zero = 0; 
    int swap_row = 0; 
    
    for (int m=l;m<=nrows;++m) {
      if (A(m,l)!=0) {
        ++non_zero;  
        swap_row = m; 
      } 
    }
    
    if (non_zero == 1) {
      if ( l!=swap_row ) {
        // swap rows
        for (int lt=1;lt<=ncols;++lt){ 
          double tmp = A(l,lt);
          A(l,lt) = A(swap_row,lt); 
          A(swap_row,lt) = tmp; 
        }
      }
    } else if (non_zero > 1) {
      // exit case 
      return l;  
    }
  }
  
  return maxcol+1; 
}


/** Compute a priori predicted estimate [z0, R0] on first sensor. Returns the predicted state vector x of the 
 *  at the position of the first sensor plane in the telescope.
 */
HepMatrix TBKalmanSRIF::ComputeBeamConstraint( HepMatrix& z0, HepMatrix& R0, ReferenceFrame& FirstSensorFrame, double mass, double mom, double charge)
{
  
  int ierr = 0; 
   
  // First, we must construct a reference uvw plane for the
  // collimator opening 
  
  ReferenceFrame CollimatorFrame;
  
  // Collimator plane is placed on z axis (beam axis), 20mm (safety margin) 
  // upstream from first sensor. 
  HepVector CollimatorPosition(3);
  CollimatorPosition[0] = 0;  
  CollimatorPosition[1] = 0;
  CollimatorPosition[2] = FirstSensorFrame.GetZPosition()-20; 
  CollimatorFrame.SetPosition(CollimatorPosition); 
  
  //cout << "Coll pos " <<   CollimatorPosition << endl; 
  
  // Collimator plane is perp. to the z axis.  
  HepMatrix CollimatorRotation;
  FillRotMatrixKarimaki(CollimatorRotation, 0,0,0);
  CollimatorFrame.SetRotation(CollimatorRotation); 
   
  //cout << "Coll rot " << CollimatorRotation << endl; 
  
  // COMPUTE REFERENCE STATE AT COLLIMATOR 
         
  HepMatrix xref_col(4,1,0);

  // COMPUTE REFERENCE STATE AT FIRST SENSOR 
  
  bool error_first = false;
  HepMatrix xref_first = TrackModel->Extrapolate(xref_col, CollimatorFrame, FirstSensorFrame, error_first);
  if (error_first) {
    streamlog_out(ERROR) << "ERR: Propagation to next surface failed. Quit fitting!"
                         << std::endl;      
  }
  
  // BEAM COVARIANCE MATRIX AT COLLIMATOR FRAME 
  
  HepSymMatrix C0(4,0);
  C0[0][0] = m_divx*m_divx;
  C0[1][1] = m_divy*m_divy; 
  C0[2][2] = m_sizex*m_sizex;
  C0[3][3] = m_sizey*m_sizey;
  C0[0][2] = m_corrx*m_divx*m_sizex;
  C0[1][3] = m_corry*m_divy*m_sizey;
  
  // COMPUTE SQUARE ROOT OF BEAM INFORMATION MATRIX
  
  C0.invert(ierr);
  R0 = diagonalize(&C0);   
   
  HepSymMatrix D0(4,0); 
  for (int i=0; i<4; ++i ) {
    D0[i][i] = std::sqrt(C0[i][i]);
  }
  
  R0 = R0*D0;  
  R0 = R0.T(); 
  
  //cout << "initial z0 " << z0 << endl;
  //cout << "intial R0 " << R0 << endl;  
  
  // MAP [z0,R0] from collimator to firsr sensor 
  // ---------------------------------------
  double theta2_col = materialeffect::GetScatterTheta2(xref_col, 0.01 , 300, mass, charge );  // 10micron kapton
  ierr = MAP_FORWARD( theta2_col, xref_col, CollimatorFrame, xref_first, FirstSensorFrame, z0, R0 );
  if (ierr != 0) {
    streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                         << std::endl;
    //return -1;
  }	
          
  //cout << "bc z0 " << z0 << endl;
  //cout << "bc R0 " << R0 << endl;  
  
  //cout << "using bc" << endl; 
  
  return xref_first;
}


} // Namespace

