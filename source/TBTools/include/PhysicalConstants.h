#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H 1

namespace depfet {

//
// System of units
//
//  Basic units:
//
//    charge        in electrons         [e]
//    distance      in milimetres        [mm]
//    energy        in electronvolts     [eV]
//    mag. field    in Tesla             [T]
//    temperature   in Kelvin            [K]
//    time          in seconds           [s]
//    voltage       in volts             [V]
//
// author: B.Schwenker, Göttingen University
//

// Elementary charge
   static const float       e = 1.;

// Charge
   static const float       C = 1/1.602176462E-19f*e;
   static const float      fC = C / 1.E15f;

// Distance
   static const float      mm = 1.;
   static const float      cm = mm * 10  ;
   static const float       m = mm * 1.E3f;
   static const float      um = mm / 1.E3f;

// Energy
   static const float      eV = 1.;
   static const float     keV = eV * 1.E3f;
   static const float     MeV = eV * 1.E6f;
   static const float     GeV = eV * 1.E9f;

// Temperature
   static const float       K = 1.;

// Time
   static const float       s = 1.;
   static const float      ms = s / 1.E3f;
   static const float      us = s / 1.E6f;
   static const float      ns = s / 1.E9f;

// Voltage
   static const float       V = 1.;

// Magnetic field
   static const float       T = 1.f*V*s/m/m;

//
//  Basic constants:
//
//    Boltzmann constant [eV/K]
//    Energy needed for creation of 1 e-h pair [eV]
//    Pi                 [1]

// Pi
   static const double  pi = 3.14159265358979323846;
   static const double piHalf = pi/2;

// Boltzmann constant in eV/K
   static const double  k = 8.617343 * 1.E-5 * eV/K;

// Energy needed for e-h pair creation
   static const double Eeh = 3.65  * eV;

// Permittivity of Si 
   static const double Perm_Si = 11.9 * 8.8542 *1.E-18 *C/V/um;

// Thermal Voltage at room temperature 
   static const double Utherm = 0.026*V;
   
// Electron mobility in intrinsic Si at room temperature 
   static const double e_mobility = 1415. * cm*cm/V/s;
 

} // Namespace

#endif // PHYSICALCONSTANTS_H
