#ifndef MaterialEffect_H
#define MaterialEffect_H 1



#include <Eigen/Core>
typedef Eigen::Matrix<double,5,1> TrackState;


namespace materialeffect {

  /** Simulate fractional energy loss due to Bremsstrahlung (Bethe Heitler theory)  
  *
  * Samples energy loss fraction z using a uniform random number rndm. The State 
  * parameter gets updated. The thickness t is given in units of the 
  * radiation length X0.
  */ 
  void SimulateBetherHeitlerEnergyLoss(TrackState& State, double t, double mass, double charge, double rndm);
  
  /** Calculate variance of the projected angular deflection due to multiple
  *  scattering of a particle for x/x0 a la Highland. Returns the scattering
  *  variance. 
  *   
  *  Projected means that the angle is not the angle in space, but rather a
  *  one-dimensional projection on one axis of a plane that is perpendicular
  *  to the particle direction before scattering.  
  */ 
  double GetScatterTheta2(const TrackState& State, double x, double x0, double mass, double charge);
  
  /** Calculate variance of the projected angular deflection due to multiple
   *  scattering of a particle for x/x0 a la Highland.
   *
   *  Projected means that the angle is not the angle in space, but rather a
   *  one-dimensional projection on one axis of a plane that is perpendicular
   *  to the particle direction before scattering. 
   */ 
  double GetScatterTheta2(double mom, double x, double x0, double mass, double charge);  
  
  /** Simulate scatter kink for single scattering theory 
   *
   * Simulate the scattering kink as a sum of independent and identically distributed single scattering 
   * events following the paper R. Frühwirth et al.  "On the quantitative modelling of core and tails of multiple 
   * scattering by Gaussian mixtures", Nucl.Instrum.Meth. (2000)
   */   
  double GetScatterKink_SC(double length, double X0, double Z, double A, double mass, double charge, double mom  );
  
  /** Scatter track at thin scatterer 
   *
   * The track is locally scattered by two scatter kink angles. The kink angles 
   * are defined as projections int the comoving frame, i.e. relative to the 
   * direction of the unscattered track. The offset through scattering is 
   * neglected (thin scatterer). The State gets updated. 
   * 
   * Note: Initially, the track state is before scattering, and gets overwritten 
   * by state after scattering.
   */ 
  void ScatterTrack(TrackState& State, double kink_u, double kink_v);
  
  /** Simulate energy loss in silicon 
   *
   * Samples an energy loss of a heavy charged particle in silicon of given 
   * thickness according to Landau theory (see PDG)
   */ 
  double SimulateEnergyLossInSilicon(TrackState& State, double thick, double mass, double charge, double lambda);
  
  /** Get most probable energy loss in silicon 
   *
   * Computes the most probable energy loss of a heavy charged particle in silicon of given 
   * thickness according to Landau theory (see PDG)
   */ 
  double GetMostProbableEnergyLossInSilicon(TrackState& State, double thick, double mass, double charge);
  
  /** Compute density effect in silicon
   *
   * Compute the density correction term for the most probable energy loss in 
   * silicon
   */ 
  double ComputeDeltaInSilicon(double gamma2, double beta2);
  
  // Highland radiation length for air [mm] (NTP: 20 deg, 1bar, PDG value)   
  static const float X0_air = 303900;

  // Average Z of air can be calculated from the fractions of the individual elements
  // According to the PDG it contains 76% nitrogen, 23% oxygen and 1% argon
  // Therefore Z=0.76*7.0+0.23*8.0+0.01*18.0
  static const float AtomicNumber_air = 7.34;

  // According to the PDG <Z/A> of air is 0.499, therefore A=Z/0.499
  static const float AtomicMass_air = 14.71; 

  // Weight to determine the modified particle momentum
  // when bremsstrahlung energy losses are significant
  // (i.e. when traversing a high Z and thick material)
  static const float Epsilon_Weightfactor = 0.3; 
 

} // Namespace


#endif
