// Local include files
#include "ReferenceFrame.h"
#include "ThreeDModel.h"

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace marlin;


namespace depfet {

// IMPLEMENT CLASS ReferenceFrame ----------------------------------------------------------------------------------------- 

// 
// Default Constructor - Identity Transformation 
// 
ReferenceFrame::ReferenceFrame() : fPosition(Eigen::Vector3d::Zero()), fRotation(Eigen::Matrix3d::Identity())
{
}

// 
// Constructur - Define Local Sensor Reference Frame 
//  	
ReferenceFrame::ReferenceFrame(Vector3d& position, Matrix3d& rotation) : fPosition(position), fRotation(rotation)
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
Vector3d ReferenceFrame::GetU()
{
  Vector3d  localU;
  localU << 1, 0, 0;
  return  fRotation.transpose() * localU;
}


Vector3d ReferenceFrame::GetV()
{
  Vector3d  localV;
  localV << 0, 1, 0;
  return  fRotation.transpose() * localV;
}

Vector3d ReferenceFrame::GetW()
{
  Vector3d localW;
  localW << 0, 0, 1;
  return  fRotation.transpose() * localW;
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
  Matrix3d deltaRot;
  FillRotMatrixKarimaki(deltaRot, dalpha, dbeta, dgamma); 
  deltaFrame.SetRotation(deltaRot);       
       
  // Shift of sensor center by dx,dy,dz in global 
  // coord. system. 
  Vector3d global_offset;
  global_offset << dx, dy, dz;

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

  Matrix3d rot = GetRotation();
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

