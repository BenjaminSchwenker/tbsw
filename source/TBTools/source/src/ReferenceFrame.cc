// Local include files
#include "ReferenceFrame.h"
#include "ThreeDModel.h"

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace CLHEP;
using namespace marlin;


namespace depfet {

// IMPLEMENT CLASS ReferenceFrame ----------------------------------------------------------------------------------------- 

// 
// Default Constructor - Identity Transformation 
// 
ReferenceFrame::ReferenceFrame() : fPosition(3, 0), fRotation(3, 3, 1)
{
}

// 
// Constructur - Define Local Sensor Reference Frame 
//  	
ReferenceFrame::ReferenceFrame(HepVector& position, HepMatrix& rotation) : fPosition(position), fRotation(rotation)
{
}

//
// Comparison operator
//
bool ReferenceFrame::operator==(const ReferenceFrame& other) const
{
  return (fPosition == other.fPosition && fRotation == other.fRotation ) ;
} 

// 
// Get direction cosines of local ref. unit vectors
//  
CLHEP::HepVector ReferenceFrame::GetU() 
{
  CLHEP::HepVector localU(3,0);
  localU[0] = 1;
  return  fRotation.T() * localU;       
}


CLHEP::HepVector ReferenceFrame::GetV() 
{
  CLHEP::HepVector localV(3,0);
  localV[1] = 1;
  return  fRotation.T() * localV;       
}

CLHEP::HepVector ReferenceFrame::GetW() 
{
  CLHEP::HepVector localW(3,0);
  localW[1] = 1;
  return  fRotation.T() * localW;  
}

ReferenceFrame ReferenceFrame::combine_karimaki(ReferenceFrame& first, ReferenceFrame& delta)
{
  ReferenceFrame combinedFrame;
  combinedFrame.SetRotation(delta.GetRotation()*first.GetRotation());
  combinedFrame.SetPosition(delta.GetPosition()+first.GetPosition());
  return combinedFrame;
}


ReferenceFrame ReferenceFrame::create_karimaki_delta(double dx, double dy, double dz, double dalpha, double dbeta, double dgamma)
{
  ReferenceFrame deltaFrame;
         
  // Small rotation by dalpha, dbeta, dgamma around
  // local sensor u-axis, (new) v-axis and (new) w- 
  // axis respectively. 
  HepMatrix deltaRot;
  FillRotMatrixKarimaki(deltaRot, dalpha, dbeta, dgamma); 
  deltaFrame.SetRotation(deltaRot);       
       
  // Shift of sensor center by dx,dy,dz in global 
  // coord. system. 
  HepVector global_offset(3);
  global_offset[0] = dx; 
  global_offset[1] = dy;      
  global_offset[2] = dz; 
  deltaFrame.SetPosition(global_offset);
  
  return deltaFrame;
} 


void ReferenceFrame::PrintHepMatrix() 
{
  streamlog_out(MESSAGE3) << " Frame Position [mm]:   "    << std::endl
                          << GetPosition()  << std::endl
                          << " Frame Rotation Matrix:   "    << std::endl 
                          << GetRotation()  << std::endl;
} 


void ReferenceFrame::PrintParams() 
{

  HepMatrix rot = GetRotation(); 
  double alpha, beta, gamma; 
  GetAnglesKarimaki(rot, alpha, beta, gamma); 
  streamlog_out(MESSAGE3) << " Frame Position [mm]:   "    << std::endl
                          << GetPosition()  << std::endl
                          << " Karimaki Alpha [rad]:   " << alpha  << std::endl
                          << " Karimaki Beta  [rad]:   " << beta   << std::endl 
                          << " Karimaki Gamma [rad]:   " << gamma  << std::endl 
                          << std::endl;
} 


} // Namespace;

