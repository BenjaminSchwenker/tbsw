#include "EVENT/LCIO.h"

namespace EVENT{
// ANSI-C++ requires initialization outside the class definition
// this has to be solved differently....
//
const std::string LCIO::LCEVENT = "LCEvent" ;
const std::string LCIO::LCCOLLECTION = "LCCollection" ;
const std::string LCIO::LCRUNHEADER = "LCRunHeader" ;
const char*  LCIO::MCPARTICLE = "MCParticle" ;
const std::string LCIO::SIMCALORIMETERHIT = "SimCalorimeterHit" ;
const std::string LCIO::RAWCALORIMETERHIT = "RawCalorimeterHit" ;
const std::string LCIO::CALORIMETERHIT = "CalorimeterHit" ;
const std::string LCIO::SIMTRACKERHIT = "SimTrackerHit" ;
const std::string LCIO::TPCHIT = "TPCHit" ;
const std::string LCIO::TRACKERRAWDATA = "TrackerRawData" ;
const std::string  LCIO::TRACKERDATA = "TrackerData" ;
const std::string LCIO::TRACKERPULSE = "TrackerPulse" ;
const std::string LCIO::TRACKERHIT = "TrackerHit" ;
const std::string LCIO::LCSTRVEC = "LCStrVec" ;
const std::string LCIO::LCFLOATVEC = "LCFloatVec" ;
const std::string LCIO::LCINTVEC = "LCIntVec" ;
const std::string LCIO::TRACK= "Track" ;
const std::string LCIO::CLUSTER = "Cluster" ;
const std::string LCIO::RECONSTRUCTEDPARTICLE = "ReconstructedParticle" ;
const std::string LCIO::LCRELATION = "LCRelation" ;
const std::string LCIO::LCGENERICOBJECT = "LCGenericObject" ;
const std::string LCIO::PARTICLEID = "ParticleID" ;
const std::string LCIO::VERTEX = "Vertex" ;

const std::string LCIO::CellIDEncoding = "CellIDEncoding" ;
}
