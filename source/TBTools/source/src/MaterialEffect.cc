// MaterialEffect implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// TBTools includes 
#include "MaterialEffect.h"

#include <cmath>
#include "TMath.h"

using namespace std;
using namespace CLHEP;


namespace materialeffect {

/** Most probable energy loss in silicon 
 *
 * Computes the most probable energy loss of a heavy charged particle in silicon of given 
 * thickness according to Landau theory (see PDG)
 */ 
double GetMostProbableEnergyLossInSilicon(double thick, double mass, double charge, double mom)
{
  double A = 28; 
  double Z = 14; 
  double rho = 0.2329; 
          
  double eMass = 0.5109989E-3;
  double exciteE = 173.0E-9;
           
  double m2 = mass*mass;
  double p2 = mom*mom;
  double e2 = m2 + p2;
          
  double beta2 = p2/e2;
  double gamma2 = e2/m2;
          
  double charge2 = charge*charge;
  
  // Energy loss spread in GeV
  double eSpread = 0.1536E-3*charge2*(Z/A)*rho*thick/beta2;
                  
  // Most probable energy loss (from integrated Bethe-Bloch equation)
  double mostProbableLoss = eSpread * ( std::log( 2.0 * eMass * beta2 * gamma2 * eSpread / (exciteE *exciteE ) )  - beta2 + 0.200 );
           
  // This one can be needed on output (but is not used internally)
  // double meanEnergyLoss = 2.*eSpread * ( std::log ( 2.*eMass*beta2*gama2 /excitE ) - beta2 );
          
  // Generate energy loss with Landau fluctuation
  return mostProbableLoss; 
} 


/** Simulate energy loss in silicon 
 *
 * Samples an energy loss of a heavy charged particle in silicon of given 
 * thickness according to Landau theory (see PDG)
*/ 
double SimulateEnergyLossInSilicon(double thick, double mass, double charge, double mom, double lambda) 
{
  double A = 28; 
  double Z = 14; 
  double rho = 0.2329; 
          
  double eMass = 0.5109989E-3;
  double exciteE = 173.0E-9;
           
  double m2 = mass*mass;
  double p2 = mom*mom;
  double e2 = m2 + p2;
          
  double beta2 = p2/e2;
  double gamma2 = e2/m2;
          
  double charge2 = charge*charge;
          
  // Energy loss spread in GeV
  double eSpread = 0.1536E-3*charge2*(Z/A)*rho*thick/beta2;
          
  // Most probable energy loss (from integrated Bethe-Bloch equation)
  double mostProbableLoss = eSpread * ( std::log( 2.0 * eMass * beta2 * gamma2 * eSpread / (exciteE *exciteE ) )  - beta2 + 0.200 );
           
  // This one can be needed on output (but is not used internally)
  // double meanEnergyLoss = 2.*eSpread * ( std::log ( 2.*eMass*beta2*gama2 /excitE ) - beta2 );
          
  // Generate energy loss with Landau fluctuation
  return mostProbableLoss + eSpread * lambda;  
} 

/** Calculate variance of the projected angular deflection due to multiple
 *  scattering of a particle for x/x0 a la Highland.
 *
 *  Projected means that the angle is not the angle in space, but rather a
 *  one-dimensional projection on one axis of a plane that is perpendicular
 *  to the particle direction before scattering. 
 */ 
double GetScatterTheta2(double x, double x0, double mass, double charge, double mom )
{
  // Sanity check - momentum mass 
  if (mass <= 0 || mom <= 0) return 0;
  
  // Sanity check - use unsigned lenght 
  if (x < 0) x*=-1;
  
  // Highland model does not like x=0 case
  if (x == 0) return 0; 
  
  // Sanity check - use air as default  
  if (x0 <= 0) x0 = 305000;    
  
  double RI = x/x0;   
  double Etot = std::sqrt(mom*mom + mass*mass);  
  double pBeta = mom*mom/Etot;  
       
  // Highland formula 
  // -------------------
  double SigTheta = 0.0136*(charge/pBeta)*std::sqrt(RI)*(1.+0.038*std::log(RI));   
  double SigTheta2 = SigTheta*SigTheta;
     
  return SigTheta2;
}  

/** Scatter track at thin scatterer 
 *
 * The track is locally scattered by two scatter kink angles. The kink angles 
 * are defined as projections int the comoving frame, i.e. relative to the 
 * direction of the unscattered track. The offset through scattering is 
 * neglected (thin scatterer). 
 * 
 * Note: Initially, the track state is before scattering, and gets overwritten 
 * by state after scattering.
 */ 
void ScatterTrack(HepMatrix& State, double kink_u, double kink_v)
{
 
  // Track direction after scatter in comoving frame
  Hep3Vector scatdir;   
  scatdir[0] = TMath::Tan(kink_u);    
  scatdir[1] = TMath::Tan(kink_v);    
  scatdir[2] = 1;  
  scatdir /= scatdir.mag(); // Change vector to unit length
  
  // Now, we construct a basis for the comoving frame
  // basis vectors are u_trk, v_trk, n_trk
   
  
  // n_trk is parallel to unscattered track direction
  Hep3Vector n_trk;    
  n_trk[0] = State[0][0];    
  n_trk[1] = State[1][0];    
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
  
  // This is the scattered dirction in 
  // detector frame :) 
  scatdir = CoRot * scatdir;
 
  // Overwritte State with new directions 
  State[0][0] = scatdir[0]/scatdir[2];    
  State[1][0] = scatdir[1]/scatdir[2];
    
  return;
  
}


}  // end namespace


