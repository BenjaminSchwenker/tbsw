#include <cstdio>
#include "gearimpl/Util.h"
#include "gear/SiPlanesLayerLayout.h"
#include "gear/LayerLayout.h"


namespace gear{

  std::ostream& operator<< (  std::ostream& s, const GearMgr& m ) {


    s << " ----------------------------------------------------- " << std::endl 
      << "              GearMgr :                                " << std::endl  
      << " ----------------------------------------------------- " << std::endl ;
    
    
    try{ 

      s << "   ----  DetectorName   ----   " << std::endl << std::endl 
	<< "         " <<  m.getDetectorName() 
	<< std::endl 
	<<  std::endl  ;

    } catch(UnknownParameterException &e){

      s << "     WARNING:  NOT FOUND  !     " << std::endl 
	<< "   please add  it to the <gear> element of your gear file: "
	<< std::endl
	<< std::endl
	<< "   <gear> " << std::endl
	<< "      <global detectorName=\"MyDetectorModle007\" /> " 
	<< std::endl
	<< "      <detectors> ... <detectors/>" << std::endl
	<< "   <gear/> " << std::endl
	<< std::endl ;  
    }

    try{ 
      s <<  m.getBField() <<  std::endl  ;
    } catch(UnknownParameterException &e){

      s << "   ----  BField   ----   " << std::endl 
	<< "     WARNING:  NOT FOUND  !     " << std::endl 
	<< "   please add  it to the <gear> element of your gear file, e.g.: "
	<< std::endl
	<< "   <gear> " << std::endl
	<< "      <BField type=\"ConstantBField\" x=\"0.0\" y=\"0.0\" z=\"4.0\"/> " 
	<< std::endl
	<< "      <detectors> ... <detectors/>" << std::endl
	<< "   <gear/> " << std::endl
	<< std::endl ;  
    }

       
    try{ 
      // SiPlanes parameters
      s <<  m.getSiPlanesParameters() <<  std::endl  ;
    } catch(UnknownParameterException &e){}  
    
 
    s << std::endl 
      << "   parameter sections (detectors) :  "
      << std::endl ;
    
    const StringVec& keys = m.getGearParameterKeys() ;
    
    for(unsigned int i=0 ; i < keys.size() ; ++i ) {
      
      s << std::endl 
	<< "   ------  name : "  <<  keys[i]  << " ------- "  << std::endl         ;

      const GearParameters& p = m.getGearParameters( keys[i] ) ;

      s <<  p <<  std::endl ;

    }

    return s ;
  }



  std::ostream& operator<< (  std::ostream& s, const   GearParameters& p ) {

    const StringVec& iKeys = p.getIntKeys() ;

    for(unsigned int i=0 ; i < iKeys.size() ; ++i ) {

      s << std::endl  << "    [int]       "  <<  iKeys[i]  <<  " \t " <<  p.getIntVal( iKeys[i] )  << std::endl ;
    }

    const StringVec& dKeys = p.getDoubleKeys() ;

    for(unsigned int i=0 ; i < dKeys.size() ; ++i ) {

      s << std::endl  << "    [double]    "  <<  dKeys[i]  <<  " \t " <<  p.getDoubleVal( dKeys[i] )  << std::endl ;
    }
    const StringVec& sKeys = p.getStringKeys() ;

    for(unsigned int i=0 ; i < sKeys.size() ; ++i ) {

      s << std::endl  << "    [string]    "  <<  sKeys[i]  <<  " \t " <<  p.getStringVal( sKeys[i] )  << std::endl ;
    }

    const StringVec& iVKeys = p.getIntVecKeys() ;

    for(unsigned int i=0 ; i < iVKeys.size() ; ++i ) {

      s << std::endl  << "    [IntVec]    "  <<  iVKeys[i]  <<  " \t "  ;

      const IntVec v = p.getIntVals( iVKeys[i] ) ;
      for(unsigned int j=0 ; j < v.size() ; ++j ) {
	s << v[j] << " , "   ;
      }

      s << std::endl ;
    }


    const StringVec& dVKeys = p.getDoubleVecKeys() ;

    for(unsigned int i=0 ; i < dVKeys.size() ; ++i ) {

      s << std::endl  << "    [DoubleVec] "  <<  dVKeys[i]  <<  " \t "  ;

      const DoubleVec v = p.getDoubleVals( dVKeys[i] ) ;
      for(unsigned int j=0 ; j < v.size() ; ++j ) {
	s << v[j] << " , "   ;
      }

      s << std::endl ;
    }

    const StringVec& sVKeys = p.getStringVecKeys() ;

    for(unsigned int i=0 ; i < sVKeys.size() ; ++i ) {

      s << std::endl  << "    [StringVec] "  <<  sVKeys[i]  <<  " \t "  ;

      const StringVec v = p.getStringVals( sVKeys[i] ) ;
      for(unsigned int j=0 ; j < v.size() ; ++j ) {
	s << v[j] << " , "   ;
      }

      s << std::endl ;
    }

    return s ;
  }

  std::ostream& operator<< (  std::ostream& s,  const BField& b ) {

    s << std::endl 
      << "   -----------   BField   ------- "  << std::endl  << std::endl        ;

    Vector3D bv = b.at(  Vector3D ( 0., 0., 0. ) ) ;

    s << "       field vector at origin :   Bx = " << bv.x() <<  " , By = " << bv.y() << " , Bz = " << bv.z()  << " [Tesla] " << std::endl ;

    s << dynamic_cast<const GearParameters&>( b )  ;
    
    return s ;

  }


  
  
  std::ostream& operator<< (  std::ostream& s,  const SiPlanesParameters& p ) {
    
    s << std::endl 
      << "   -----------   SiPlanesParameters  ------- "  << std::endl         ;
    
    s << dynamic_cast<const GearParameters&>( p )  ;
    
    const SiPlanesLayerLayout & l = p.getSiPlanesLayerLayout() ;


    s <<  std::endl << " Setup ID : " << p.getSiPlanesID() << std::endl;

    int type = p.getSiPlanesType() ;
    
    s << " Telescope type : " ;
    std::string strType ;

    switch( type ) {
      
    case SiPlanesParameters::TelescopeWithDUT : 
      strType = "TelescopeWithDUT" ;
      s << " TelescopeWithDUT" << std::endl ;       
      break ;
    case SiPlanesParameters::TelescopeWithoutDUT : 
      strType = "TelescopeWithoutDUT" ;
      s << " TelescopeWithoutDUT " << std::endl ;      
      break ;
      
    default :  
      
      s << " unknown " << std::endl ; 
    }

    s << " Number of telescope planes : " << p.getSiPlanesNumber() << std::endl;
    
    
    // layers
    
									      
    s << " Number of layers : " << l.getNLayers() << std::endl << std::endl ;

    s <<  " layer parameters : "  << std::endl ;

    char buffer[1024] ;
    
    std::sprintf(buffer,"  |-------------------------------------------------------------------------------------------------------------------------------------------------|\n") ;
    s << buffer ;

    std::sprintf(buffer,"  |              ladder:                        |                    sensitive part:                                                                |\n") ;
    s << buffer ;

    std::sprintf(buffer,"  | ID | pozX| pozY|  pozZ | sizeX| sizeY| Thick| ID| pozX| pozY|  pozZ |sizeX|sizeY| Thick|NpixX|NpixY|PitchX|PitchY| Resol| Rot1| Rot2| Rot3| Rot4| \n") ;

    s << buffer ;

    std::sprintf(buffer,"  |-------------------------------------------------------------------------------------------------------------------------------------------------|\n") ;
    s << buffer ;

    for( int i = 0 ; i < l.getNLayers() ; i++ ) {

      char buffer1[1024] ;
      std::sprintf(buffer1,"  |%4d|%5.2f|%5.2f|%7.2f|%6.2f|%6.2f|%6.3f|%3d|%5.2f|%5.2f|%7.2f|%5.2f|%5.2f|%6.3f| %4d| %4d|%6.2f|%6.2f|%6.4f| %4.2f| %4.2f| %4.2f| %4.2f|\n"
	      , l.getID(i) 
	      , l.getLayerPositionX(i) 
	      , l.getLayerPositionY(i) 
	      , l.getLayerPositionZ(i)
	      , l.getLayerSizeX(i) 
	      , l.getLayerSizeY(i) 
	      , l.getLayerThickness(i)
  	      , l.getSensitiveID(i) 
 	      , l.getSensitivePositionX(i) 
	      , l.getSensitivePositionY(i) 
	      , l.getSensitivePositionZ(i)
	      , l.getSensitiveSizeX(i) 
	      , l.getSensitiveSizeY(i) 
	      , l.getSensitiveThickness(i) 
	      , l.getSensitiveNpixelX(i) 
	      , l.getSensitiveNpixelY(i) 
	      , l.getSensitivePitchX(i) 
	      , l.getSensitivePitchY(i) 
	      , l.getSensitiveResolutionX(i)
	      , l.getSensitiveRotation1(i) 
	      , l.getSensitiveRotation2(i) 
	      , l.getSensitiveRotation3(i) 
	      , l.getSensitiveRotation4(i) ) ;
      
      s << buffer1 ;

    }

    std::sprintf(buffer,"  |-------------------------------------------------------------------------------------------------------------------------------------------------|\n") ;
    s << buffer ;

    // DUT

    if (strType == "TelescopeWithDUT"){

      s <<  " DUT parameters : "  << std::endl ;
      
      s << buffer ;
      
      std::sprintf(buffer,"  |              ladder:                        |                    sensitive part:                                                         |\n") ;
      s << buffer ;
      
      std::sprintf(buffer,"  |------------------------------------------------------------------------------------------------------------------------------------------|\n") ;
      s << buffer ;
      
      std::sprintf(buffer,"  | ID | pozX| pozY|  pozZ | sizeX| sizeY| Thick| ID| pozX| pozY|  pozZ |sizeX|sizeY| Thick|NpixX|NpixY|PitchX|PitchY| Rot1| Rot2| Rot3| Rot4| \n") ;
      s << buffer ;
      
      std::sprintf(buffer,"  |------------------------------------------------------------------------------------------------------------------------------------------|\n") ;
      s << buffer ;
      
      char buffer1[1024] ;
      std::sprintf(buffer1,"  |%4d|%5.2f|%5.2f|%7.2f|%6.2f|%6.2f|%6.3f|%3d|%5.2f|%5.2f|%7.2f|%5.2f|%5.2f|%6.3f| %4d| %4d|%6.2f|%6.2f| %4.2f| %4.2f| %4.2f| %4.2f|\n"
	      , l.getDUTID() 
	      , l.getDUTPositionX() 
	      , l.getDUTPositionY() 
	      , l.getDUTPositionZ()
	      , l.getDUTSizeX() 
	      , l.getDUTSizeY() 
	      , l.getDUTThickness() 
	      , l.getDUTSensitiveID() 
	      , l.getDUTSensitivePositionX() 
	      , l.getDUTSensitivePositionY() 
	      , l.getDUTSensitivePositionZ()
	      , l.getDUTSensitiveSizeX()
	      , l.getDUTSensitiveSizeY() 
	      , l.getDUTSensitiveThickness() 
	      , l.getDUTSensitiveNpixelX() 
	      , l.getDUTSensitiveNpixelY() 
	      , l.getDUTSensitivePitchX() 
	      , l.getDUTSensitivePitchY() 
	      , l.getDUTSensitiveRotation1() 
	      , l.getDUTSensitiveRotation2() 
	      , l.getDUTSensitiveRotation3() 
	      , l.getDUTSensitiveRotation4() ) ;
      
      s << buffer1 ;
      
      std::sprintf(buffer,"  |------------------------------------------------------------------------------------------------------------------------------------------|\n") ;
      s << buffer ;
      
    }
    return s ;
    
    
  }

}
