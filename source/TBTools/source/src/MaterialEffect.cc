// MaterialEffect implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


// TBTools includes 
#include "MaterialEffect.h"

#include <cmath>
#include "TMath.h"
#include <Math/DistFunc.h>

using namespace std;
using namespace CLHEP;


namespace materialeffect {

/** Simulate fractional energy loss due to Bremsstrahlung (Bethe Heitler theory)  
 *
 * Samples energy loss fraction z using a uniform random number rndm. The State 
 * parameter gets updated. The thickness t is given in units of the 
 * radiation length X0.
 */ 
void SimulateBetherHeitlerEnergyLoss(HepMatrix& State, double t, double mass, double charge, double rndm)     
{
  if ( t > 0  ) { 
    // Comput energy loss fraction
    double c = t/TMath::Log(2);
    double u = ROOT::Math::gamma_quantile(rndm,c,1.0);
    double z = TMath::Exp(-u);
    
    // Compute initial energy 
    double mom = charge/State[4][0];
    double m2 = mass*mass;
    double p2 = mom*mom;
    double e2 = m2 + p2;
    
    // Compute momentum after energy loss 
    e2 *= z*z; 
    p2 = e2 - m2;
    if (p2>0) mom = TMath::Sqrt(p2);
    else mom = 0.001*mass;
    
    // Finally, update the state
    State[4][0] = charge/mom;  
  } 
}

/** Compute density effect in silicon
 *
 * Compute the density correction term for the most probable energy loss in 
 * silicon
 */ 
double ComputeDeltaInSilicon(double gamma2, double beta2)
{

  double S0 = 0.201;
  double S1 = 2.872;
  double hv = 31.06E-9;
  double md = 3.255;
  double a = 0.149;
  double delta0 = 0.14;
  double exciteE = 170.0E-9;
  
  double betagamma = sqrt(gamma2*beta2);
  double A = a*std::pow( log( pow(10.0,S1) / betagamma )/log(10.0) ,md);
  double C = -2.0*log(exciteE/hv) -1 ; 
  
  if (betagamma < std::pow(10.0,S0)) {
    return delta0*std::pow(betagamma/std::pow(10.0,S0),2);
  } else if ( betagamma >  std::pow(10.0,S1) ) { 
    return  2*log(sqrt(beta2*gamma2)) + C;
  } 
  return  2*log(sqrt(beta2*gamma2)) + C + A;
}




/** Most probable energy loss in silicon 
 *
 * Computes the most probable energy loss of a heavy charged particle in silicon of given 
 * thickness according to Landau theory (see PDG)
 */ 
double GetMostProbableEnergyLossInSilicon(HepMatrix& State, double thick, double mass, double charge)
{
  double A = 28; 
  double Z = 14; 
  double rho = 0.2329; 
          
  double eMass = 0.5109989E-3;
  double exciteE = 170.0E-9;
         
  double mom = charge/State[4][0];   
  double m2 = mass*mass;
  double p2 = mom*mom;
  double e2 = m2 + p2;
          
  double beta2 = p2/e2;
  double gamma2 = e2/m2;
          
  double charge2 = charge*charge;

  // Density correction 
  double delta =  ComputeDeltaInSilicon(gamma2, beta2); 
  
  // Energy loss spread in GeV
  double eSpread = 0.1536E-3*charge2*(Z/A)*rho*thick/beta2;
                  
  // Most probable energy loss (from integrated Bethe-Bloch equation)
  double mostProbableLoss = eSpread * ( std::log( 2.0 * eMass * beta2 * gamma2 * eSpread / (exciteE *exciteE ) )  - beta2 + 0.200 - delta*0.5);
                    
  // Generate energy loss with Landau fluctuation
  return mostProbableLoss; 
} 


/** Simulate energy loss in silicon 
 *
 * Samples an energy loss of a heavy charged particle in silicon of given 
 * thickness according to Landau theory (see PDG)
*/ 
double SimulateEnergyLossInSilicon(HepMatrix& State, double thick, double mass, double charge, double lambda) 
{
  double A = 28; 
  double Z = 14; 
  double rho = 0.2329; 
          
  double eMass = 0.5109989E-3;
  double exciteE = 170.0E-9;
           
  double mom = charge/State[4][0]; 
  double m2 = mass*mass;
  double p2 = mom*mom;
  double e2 = m2 + p2;
          
  double beta2 = p2/e2;
  double gamma2 = e2/m2;
          
  double charge2 = charge*charge;
          
  // Density correction 
  double delta =  ComputeDeltaInSilicon(gamma2, beta2); 
  
  // Energy loss spread in GeV
  double eSpread = 0.1536E-3*charge2*(Z/A)*rho*thick/beta2;
          
  // Most probable energy loss (from integrated Bethe-Bloch equation)
  double mostProbableLoss = eSpread * ( std::log( 2.0 * eMass * beta2 * gamma2 * eSpread / (exciteE *exciteE ) )  - beta2 + 0.200 - delta*0.5);
                   
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
double GetScatterTheta2(HepMatrix& State, double x, double x0, double mass, double charge)  
{
  double mom = charge/State[4][0]; 
  
  // Sanity check -  mass 
  if (mass <= 0) return 0;
  
  // Sanity check - momentum 
  if (mom <= 0) return 0;

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


/** Calculate variance of the projected angular deflection due to multiple
 *  scattering of a particle for x/x0 a la Highland.
 *
 *  Projected means that the angle is not the angle in space, but rather a
 *  one-dimensional projection on one axis of a plane that is perpendicular
 *  to the particle direction before scattering. 
 */ 
double GetScatterTheta2(double mom, double x, double x0, double mass, double charge)  
{
  // Sanity check -  mass 
  if (mass <= 0) return 0;
  
  // Sanity check - momentum 
  if (mom <= 0) return 0;

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


