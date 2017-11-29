// TBKalmanMSC implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// Marlin includes 
#include <marlin/Global.h>
#include <streamlog/streamlog.h>
  
// DEPFETTrackTool includes
#include "TBKalmanMSC.h"
#include "MaterialEffect.h"
#include "StraightLineTrackModel.h"
#include "HelixTrackModel.h"

// C++ STL includes 
#include <cmath>
#include <limits>

// Namespaces
using namespace CLHEP;
using namespace marlin;
using namespace std; 

namespace depfet {

/** Constructor 
 */
TBKalmanMSC::TBKalmanMSC(TBDetector& detector)
{                   
  
  ndim = 5; // dimension of state vector

  // Project hit coord out of track state
  H = HepMatrix(2,ndim,0);
  H[0][2] = 1.;
  H[1][3] = 1.;

  mom = 0; 
  mass = 0;
  charge = 0;

  HepVector field(3,0);
  field[0] = detector.GetBx();
  field[1] = detector.GetBy();
  field[2] = detector.GetBz(); 
  
  if ( field.norm() == 0 ) {
    TrackModel = new StraightLineTrackModel();     
  } else {
    TrackModel = new HelixTrackModel(field); 
  }
}

TBKalmanMSC::~TBKalmanMSC()
{ 
  delete TrackModel; 
}


/** Set number of degrees of freedom in trk
 */
void TBKalmanMSC::SetNdof(TBTrack& trk)
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


/** Performs track fitting. Returns track chisqu. 
 */
bool TBKalmanMSC::ProcessTrack(TBTrack& trk, int dir, bool biased)
{
  
  // Get particle hypothesis 
  mom = trk.GetMomentum(); 
  mass = trk.GetMass();
  charge = trk.GetCharge();
   
  // Get all track elements (TEs) in track
  std::vector<TBTrackElement>& TEVec = trk.GetTEs();
  int nTE = (int) TEVec.size();
   
  int IFSURF = 0;      // First TE 
  int ILSURF = nTE-1;  // Last TE  
  
  // The bare minimum for a meaningfull track state 
  // estimate are two hits. For three or more hits, 
  // the track chisqu > 0 and may be used for tracking
  // cuts. 
     
  // Iterations of double Kalman filter    
  double chisqu = -1;    
  trk.SetChiSqu(chisqu);
  SetNdof(trk);     
  
  if ( trk.GetNumHits() < 2 ) {  
    trk.SetChiSqu(-1);  
    return true;
  }   

  for(int iter=0; iter<2; iter++) { 
    
    // Fit linearization
    // 
    // Compute intersections of the reference trajectory
    // with sub detectors.  
     
    // Get paramters of reference trajectory   
    HepMatrix RTrkState = trk.GetReferenceState().GetPars();      
    ReferenceFrame RTrkSurf = trk.GetTE( trk.GetReferenceState().GetPlane() ).GetDet().GetNominal();
     
    REFTrack RefStateVec(nTE);
    
    for(int iTE=0; iTE<nTE; ++iTE) {
       
      // Next surface along beam line  
      const ReferenceFrame& Next_Surf = TEVec[iTE].GetDet().GetNominal();
           
      // Get reference state on Next_Surf
      bool error = false;
      RefStateVec[iTE] = TrackModel->Extrapolate(RTrkState, RTrkSurf, Next_Surf, error);
      if (error) {
        streamlog_out(ERROR) << "ERR: Msc Propagation to next detector failed. Quit fitting!"
                             << std::endl;
        return error;       
      }
      
    } 
    
    // Initialize Kalman filter 
    
    // This is the initial deviation between track state and 
    // reference track   
    HepMatrix x0(ndim,1,0);
     
    // This is the covariance matrix of the above 
    // deviation vector 
    HepSymMatrix C0(ndim,0);
    
    for (int i = 0; i < 2; i++) {
      // There is a tradeoff between choosing large values (which makes the
      // tracking numerically unstable) and choosing small values (which biases
      // the track from the seed). I found 1E4 working for this model, but this
      // might need to get tuned.
      C0[i][i] = 1E-2;
      C0[i+2][i+2] = 1E1; 
    } 
    //C0[4][4] = 100;  
    C0[4][4] = 1;  
    
    // Store for filter outputs 
    KalFilterVecMSC KalmanFilterData(nTE);
     
    if ( dir < 0 ) { 
      // Backward KALMAN filter -----------------------------------------------------
      KalmanFilterData[ILSURF].Pr_x = x0;
      KalmanFilterData[ILSURF].Pr_C = C0;
      chisqu = KalmanFilter(trk, -1, ILSURF, IFSURF, KalmanFilterData, RefStateVec);
    } else {
      // Forward KALMAN filter -----------------------------------------------------
      KalmanFilterData[IFSURF].Pr_x = x0;
      KalmanFilterData[IFSURF].Pr_C = C0;
      chisqu = KalmanFilter(trk, +1, IFSURF, ILSURF, KalmanFilterData, RefStateVec);    
    }
    
    // A very basic consistency test 
    if ( std::isnan(chisqu) ||  chisqu < 0 )  {
      trk.SetChiSqu(-1);
      return true;
    }
      
    // Store results 
    
    for(int iTE=0; iTE<nTE; ++iTE) {
       
      TBTrackElement& TE = TEVec[iTE];
      TE.SetChiSqu(KalmanFilterData[iTE].Chi2Pred); 
       
      if (!biased) { 
        TE.GetState().GetPars() = KalmanFilterData[iTE].Pr_x;  
        TE.GetState().GetCov() =  KalmanFilterData[iTE].Pr_C;  
      } else {
        TE.GetState().GetPars() = KalmanFilterData[iTE].Up_x;  
        TE.GetState().GetCov() = KalmanFilterData[iTE].Up_C; 
      }
       
      // Do not forget to add reference state         
      TE.GetState().GetPars() += RefStateVec[iTE];
    }

    
    // Relinearize the track fit for next iteration around  
    // smoothed state at IFSURF
         
    // Use full information from fit for reference state 
    int ipl; 
    if (dir<0) {
      ipl = ILSURF;     
    } else {
      ipl = IFSURF;     
    } 
     
    trk.GetReferenceState().Pars = trk.GetTE( ipl ).GetState().GetPars();  
    trk.GetReferenceState().SetPlane(ipl);  
    trk.SetMomentum(std::abs(charge/trk.GetTE( ipl ).GetState().GetPars()[4][0]));
    
    
  }
  
  trk.SetChiSqu(chisqu);
     
  // Everything is ok
  return false;
}

CLHEP::HepMatrix TBKalmanMSC::GetScatterKinks(Det& DetUnit, TBTrackState& InState, TBTrackState& OutState)
{
  double slopes[4];
  slopes[0]=InState.GetPars()[0][0];
  slopes[1]=InState.GetPars()[1][0];
  slopes[2]=OutState.GetPars()[0][0];
  slopes[3]=OutState.GetPars()[1][0]; 
  return slopestotheta(slopes);  
} 

CLHEP::HepSymMatrix TBKalmanMSC::GetScatterKinkCov(Det& DetUnit, TBTrackState& InState, TBTrackState& OutState)
{

  //define aid varibales
  double Slope[4];
  double Cov[4];
  
  //weights used for the unscented transform
  double meanweight[9];
  double covweight[9];
  
  //Sigmapoints used in the unscented transform
  double Sigmapoint[9][4];
  
  //Theta values of the Sigmapoints
  HepMatrix Theta(2,9,0);

  //reseting the mean of the unscented transform
  double mean1=0;
  double mean2=0;

  //reseting the covariance matrix diagonal elements
  double cov12=0;
  double cov1=0;
  double cov2=0;

  //defining the weight parameters
  double alpha=TMath::Sqrt(3);   //parameter used in the weight of the covariance matrix
  double beta=2;                 //parameter used in the weight of the covariance matrix
  double lambda=1;               //scaling parameter
  double n=4;                    //dimensionality
  
  //Get the Slopes and Covs from the In and OutState
  HepMatrix State_in_pars = InState.GetPars();
  HepMatrix State_out_pars = OutState.GetPars();
  HepMatrix State_in_covs = InState.GetCov();
  HepMatrix State_out_covs = OutState.GetCov();

  Slope[0]=State_in_pars[0][0];
  Slope[1]=State_in_pars[1][0];
  Slope[2]=State_out_pars[0][0];
  Slope[3]=State_out_pars[1][0];
  
  Cov[0]=State_in_covs[0][0];
  Cov[1]=State_in_covs[1][1];
  Cov[2]=State_out_covs[0][0];
  Cov[3]=State_out_covs[1][1];
  
  //Compute first Sigmapoint
  for(int i=0;i<4;i++) {
    Sigmapoint[0][i]=Slope[i];
  }
  
  //Compute second Sigmapoint
  Sigmapoint[1][0]=Slope[0]+TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[0]);
  for(int i=1;i<4;i++) {
    Sigmapoint[1][i]=Slope[i];
  }
  
  //Compute third Sigmapoint
  Sigmapoint[2][0]=Slope[0];
  Sigmapoint[2][1]=Slope[1]+TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[1]);
  for(int i=2;i<4;i++) {
    Sigmapoint[2][i]=Slope[i];
  }
	 
  //Compute fourth Sigmapoint
  for(int i=0;i<2;i++) {
    Sigmapoint[3][i]=Slope[i];
  }  
  Sigmapoint[3][2]=Slope[2]+TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[2]);
  Sigmapoint[3][3]=Slope[3];

  //Compute fifth Sigmapoint
  for(int i=0;i<3;i++) {
    Sigmapoint[4][i]=Slope[i];
  }
  Sigmapoint[4][3]=Slope[3]+TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[3]);
  
  //Compute sixth Sigmapoint
  Sigmapoint[5][0]=Slope[0]-TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[0]);
  for(int i=1;i<4;i++){
    Sigmapoint[5][i]=Slope[i];
  }

  //Compute seventh Sigmapoint
  Sigmapoint[6][0]=Slope[0];
  Sigmapoint[6][1]=Slope[1]-TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[1]);
  for(int i=2;i<4;i++) {
    Sigmapoint[6][i]=Slope[i];
  }

  //Compute eighth Sigmapoint
  for(int i=0;i<2;i++){
    Sigmapoint[7][i]=Slope[i];
  }
  Sigmapoint[7][2]=Slope[2]-TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[2]);
  Sigmapoint[7][3]=Slope[3];

  //Compute ninth Sigmapoint
  for(int i=0;i<3;i++) {
    Sigmapoint[8][i]=Slope[i];
  }
  Sigmapoint[8][3]=Slope[3]-TMath::Sqrt(n+lambda)*TMath::Sqrt(Cov[3]);
  
  //computing the weighths
  
  //first the mean weigths
  meanweight[0]=lambda/(n+lambda); 
  for(int i =1;i<9;i++)  {
    meanweight[i]=lambda/(2*(n+lambda));
  }

  //now the covariance weights
  covweight[0]=meanweight[0]+(1.0-alpha*alpha+beta); 
  for(int i=1;i<9;i++){
    covweight[i]=lambda/(2*(n+lambda));
  }
  
  //Compute the thetas for each sigma point
  for(int i=0;i<9;i++)   {
  
    HepMatrix helptheta(2,1,0);
    double help[4];
    for(int j=0;j<4;j++){
      help[j]=Sigmapoint[i][j];
    }
      
    //for the thetas
    //Fehler liegt hier ich muss für jeden slope aus den sigmapoints die thetas bestimmen
    helptheta=slopestotheta(help);
    Theta[0][i]=helptheta[0][0];
    Theta[1][i]=helptheta[1][0];	
  }
  
  //Computing the mean of the unscented transform
  for(int i=0;i<9;i++) {
    mean1=mean1+meanweight[i]*Theta[0][i];
    mean2=mean2+meanweight[i]*Theta[1][i];
  }

  //Computing the diagonal elements of the covariance matrix of the unscented transform
  for(int i=0;i<9;i++){
    cov1=cov1+covweight[i]*(Theta[0][i]-mean1)*(Theta[0][i]-mean1);
    cov2=cov2+covweight[i]*(Theta[1][i]-mean2)*(Theta[1][i]-mean2);
    cov12=cov12+covweight[i]*(Theta[0][i]-mean1)*(Theta[1][i]-mean2);
  }
  HepSymMatrix Covariance(2,0);
  Covariance[0][0]=cov1;
  Covariance[1][0]=cov12;
  Covariance[1][1]=cov2;

  return Covariance;
}





//Function for getting the thetas from the slopes
CLHEP::HepMatrix TBKalmanMSC::slopestotheta(double slopes[4]){
     
  //Scatter direction in the detector frame
  Hep3Vector scatdir;
  scatdir[0] = slopes[2];
  scatdir[1] = slopes[3];
  scatdir[2] = 1;
  scatdir /= scatdir.mag(); //Change vector to unit length
  
  // Now, we construct a basis for the comoving frame
  // basis vectors are u_trk, v_trk, n_trk
   
  // n_trk is parallel to unscattered track direction
  Hep3Vector n_trk;  
  n_trk[0] = slopes[0];    
  n_trk[1] = slopes[1];    
  n_trk[2] = 1;  
  n_trk /= n_trk.mag(); // Change vector to unit length
  
  // v_trk is orthogonal to track dir and detector u axis  
  Hep3Vector u_hat(1,0,0);
  Hep3Vector v_trk = n_trk.cross(u_hat);
  v_trk /= v_trk.mag();
  
  // u_trk completes rigth handed system  
  Hep3Vector u_trk = v_trk.cross(n_trk);   
  
  // Now, we construct rotation matrix from comoving 
  // frame to detector frame 
  HepRotation CoRot; 
  CoRot.rotateAxes(u_trk, v_trk, n_trk); 
  
  //Determine the inverese matrix for the transformation from the detector to the comoving frame
  HepRotation CoRot_inv;
  CoRot_inv = CoRot.inverse();
  
  //Get the old scatter direction in the comoving frame again by using the inverse matrix
  Hep3Vector xvec_a;
  xvec_a = CoRot_inv * scatdir;
  
  //compute the slopes from the scattering direction in the comoving frame
  double Slopex = xvec_a[0]/xvec_a[2];
  double Slopey = xvec_a[1]/xvec_a[2];
  
  //compute the angles by applying ArcusTangens on the slopes
  HepMatrix theta(2,1,0);
  theta[0][0] = TMath::ATan(Slopex);
  theta[1][0] = TMath::ATan(Slopey);
     
  return theta;
}


double TBKalmanMSC::KalmanFilter(TBTrack& trk, int idir, int istart, int istop, KalFilterVecMSC& result, REFTrack& RStateVec)
{ 
      
  // To start ierr is set to 0 (= OK)
  int ierr = 0; 
  // Final chisqu of track fit
  double chisqu = 0;
       
  

  // Process vector of track elements
  std::vector<TBTrackElement>& TEVec = trk.GetTEs();  
  
  // Note that the Kalman filter for tracking is linearized 
  // around the reference trajectory.     
  // It means that the Kalman filter is actually applied to 
  // the deviations from the reference trajectory. 
  
  for(int is=istart; is!=istop+idir; is+=idir) {
    
    //cout << "isurface "<< is << endl;
    
    // Process next track element 
    TBTrackElement& te = TEVec[is];
     
    // Get surface parameters           
    ReferenceFrame& Surf = te.GetDet().GetNominal();  
    
    // Get current reference state xref
    HepMatrix& xref = RStateVec[is];    
    
    
    // Copy predicted [x0,C0]
    HepMatrix x0 = result[is].Pr_x;  
    HepSymMatrix C0 = result[is].Pr_C;
     
    //cout << "Predicted x0 " << x0 << endl;
    //cout << "Predicted C0 " << C0 << endl;
     
    // If measurement, add information and update chi2   
    double predchi2 = 0;  
     
    if (te.HasHit()) {
             
      // Measured hit coordinates, 2x1 matrix 
      HepMatrix m = te.GetHit().GetCoord();
      
      // Covariance for hit coordinates, 2x2 matrix 
      HepSymMatrix V = te.GetHit().GetCov();
      
      // Measured deviation from the 
      // reference trajectory
      m -= H*xref;
               
      // Weigth matrix of measurment 
      HepSymMatrix W = (V + C0.similarity(H)).inverse(ierr);
      if (ierr != 0) {
        streamlog_out(ERROR) << "ERR: Matrix inversion failed. Quit fitting!"
                             << std::endl;
        return -1;
      }	
       
      // This is the predicted residual 
      HepMatrix r = m - H*x0; 
      
      // This is the predicted chi2 
      HepMatrix chi2mat = r.T()*W*r;
      predchi2 = chi2mat[0][0];   
      
      // Kalman gain matrix K 
      HepMatrix K = C0 * H.T() * W; 
       
      // This is the filtered deviation from the 
      // reference trajectory.  
      x0 += K * r;
      C0 -= (W.similarityT(H)).similarity(C0);
    }  
      
    
    
    // Store results and update chisqu
    chisqu += predchi2;
    result[is].Chi2Pred = predchi2; 
    result[is].Up_x = x0;  
    result[is].Up_C = C0;
    
    // -----------------------------    
    // Propagate to next surface
    // -----------------------------
    
    // The linearized mapping from state at k to state at k+1:
    // 
    // x_k+1 = J(k+1,k) * x_k + G_k * w_k                    (*)
    // 
    //  J(j,k) : Transport matrix from k to j -> linearized around ref state  
    //  w_k      : Vector of all scatter angles for scatterers between k->k+1
    //  Q_k      : MSC covariance matrix for w_k -> w_k are zero mean variates 
    //  G_k      : Scatter gain matrix -> linearized around ref state   
    // 
    // The first part in equation (*) is just a straight line extrapolation, 
    // but the second part includes the influence of random scatterings w_k
    // along the way. In particular, w_k includes the scatterings at sensor 
    // k itself. From a physics perspective, x_k estimates 'in' states, i.e. 
    // the track parameters before the scattering in sensor k happens.
    // 
    // Note that eq. (*) defines a forward mapping with process noise. In 
    // HEP tracking, transport matrix is typically invertible. The backward 
    // mapping looks 
    // 
    // x_k = J(k+1,k)^-1 * x_k+1 - J(k+1,k)^-1 * G_k * w_k   (**) 
    // 
    // To study multiple scattering processes, it makes sense to develop estimators 
    // for the out states as well. To do so, we revers the time direction in eq. (*) 
    // and find the transport equation for out states 
    // 
    // xo_k-1 = J(k-1,k) * xo_k + G_k * w_k                  (***)
    // 
    // Now, we have estimator for the track direction before and after scattering for 
    // all inner sensors. Please note that it is ill defined to calculate a average 
    // between 'in' and 'out' states. They life in different spaces. The estimation of 
    // scatter angles can be done by computing a difference between in and out states.      
    
    // Plane number of next surface  
    int inext = is+idir;
     
    if (inext!=istop+idir)  {
      
      // Next surface along filter direction 
      TBTrackElement& nte = TEVec[inext]; 
      
      // Parameters for next surface          
      ReferenceFrame& nSurf = nte.GetDet().GetNominal();
      
      // Transport matrix
      HepMatrix J(ndim, ndim);
      ierr = TrackModel->TrackJacobian( xref, Surf, nSurf, J);
      if (ierr != 0) {
        streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                             << std::endl;
        return -1;
      }	
      
      // Time update of state    
      HepMatrix x1 = J*x0; 
      
      // Time update of covariance 
      // 
      // This fitter aims to take into account scatter in 
      // air between sensors. Therefor, we add a virtual 
      // air surface between the detectors which scatters 
      // the track with a material budget equal to the 
      // extrapolation step length between the sensors.  
              
      // Get signed flight length in air between detectors 
      double length = TrackModel->GetSignedStepLength(xref, Surf, nSurf);
      
      // And extrapolate reference track to air scatter surface
      HepMatrix xref_air = xref; 
      ReferenceFrame Surf_air = Surf; 
      TrackModel->Extrapolate(xref_air, Surf_air, length/2);
            
      // Extrap covariance to air surface
      // -------------------------------- 
      HepMatrix J_det2air(ndim, ndim);
      ierr = TrackModel->TrackJacobian( xref, Surf, Surf_air, J_det2air);
      if (ierr != 0) {
        streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                             << std::endl;
        return -1;
      }	
      
      HepSymMatrix C1 = C0.similarity(J_det2air);
      
      // Add scatter noise from last detector
      // -----------------------------------
      
      // Local Scatter gain matrix  
      HepMatrix Gl_det(ndim, 2);    
      TrackModel->GetScatterGain(xref, Gl_det);
          
      // General Scatter gain matrix 
      HepMatrix G_det = J_det2air*Gl_det;   
        
      // Variance of projected scatter angles   
      HepSymMatrix Q_det(2, 1);    

      double l0 = te.GetDet().GetTrackLength(xref[2][0], xref[3][0], xref[0][0], xref[1][0]);
      double X0 = te.GetDet().GetRadLength(xref[2][0],xref[3][0]);    
      double theta2_det = materialeffect::GetScatterTheta2(xref, l0, X0, mass, charge);    
      Q_det *= theta2_det;
            
      C1 += Q_det.similarity(G_det);    
       
      // Extrap covariance to next detector
      // ----------------------------------    
      HepMatrix J_air2det(ndim, ndim);
      ierr = TrackModel->TrackJacobian( xref_air, Surf_air, nSurf, J_air2det);
      if (ierr != 0) {
        streamlog_out(ERROR) << "ERR: Problem with track extrapolation. Quit fitting!"
                             << std::endl;
        return -1;
      }	
      
      C1 = C1.similarity(J_air2det);
      
      // Add scatter noise from air 
      //---------------------------   
      
      // Local Scatter gain matrix  
      HepMatrix Gl_air(ndim, 2);    
      TrackModel->GetScatterGain(xref_air, Gl_air);
          
      // General Scatter gain matrix 
      HepMatrix G_air = J_air2det*Gl_air;   
        
      // Variance of projected scatter angles   
      HepSymMatrix Q_air(2, 1);     
      double theta2_air = materialeffect::GetScatterTheta2(xref_air, length, materialeffect::X0_air, mass, charge );   
      Q_air *= theta2_air;
            
      C1 += Q_air.similarity(G_air);    
       
      
      // Store prediction results        
      result[inext].Chi2Pred = 0; 
      result[inext].Pr_x = x1;  
      result[inext].Pr_C = C1;
      
    } 
           
  }
      
  // Successfull return  
  return chisqu;  
}


} // Namespace;

