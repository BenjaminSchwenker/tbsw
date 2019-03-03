// Det implementation file 
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

// TBTools includes 
#include "Det.h"

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
using namespace marlin;

namespace depfet {



Det::Det(const std::string& typeName, int sensorID, int planeNumber) :  
  m_typeName(typeName),
  m_sensorID(sensorID),
  m_planeNumber(planeNumber), 
{
  //register detector in map
  //DetMgr::instance()->registerDet( this ) ;
}

Det::~Det() {}


/**  Print methods
 */ 
void Det::Print() const
{
  streamlog_out(MESSAGE3) << std::endl
                          << "  Type     :    " << m_typeName    << std::endl
                          << "  Plane no.:    " << m_planeNumber << std::endl  
                          << "  Sensor ID:    " << m_sensorID    << std::endl;  
}
 


} // Namespace;

