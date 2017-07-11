// ClusterShapeToAsciiPrinter Processor  
// 
// See ClusterShapeToAsciiPrinter.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>


#include "ClusterShapeToAsciiPrinter.h"

// TBTools includes
#include "TBTrack.h"
#include "TrackInputProvider.h"
#include "GenericTrackFitter.h"
#include "PixelCluster.h"

// Include basic C
#include <iostream>
#include <limits>
#include <iomanip>

// Include LCIO classes
#include "lcio.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <Exceptions.h>

// Used namespaces
using namespace std; 
using namespace CLHEP; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

  //
  // Instantiate this object
  //
  ClusterShapeToAsciiPrinter aClusterShapeToAsciiPrinter ;
  
  //
  // Constructor
  //
  ClusterShapeToAsciiPrinter::ClusterShapeToAsciiPrinter() : Processor("ClusterShapeToAsciiPrinter")
  {
    
    // Processor description
    _description = "ClusterShapeToAsciiPrinter: Print cluster shapes to ascii file" ;
    
    //
    // Input collections 
    
    registerInputCollection(LCIO::TRACK,"InputTrackCollectionName",
                            "Track input collection",
                            _inputTrackCollectionName,std::string("tracks"));
    
    registerProcessorParameter( "OuputFileName",
                                "File name containing cluster shapes",
                                _asciiFileName, std::string("ClusterShapes.txt"));  
      
    registerProcessorParameter ("AlignmentDBFileName",
                                "This is the name of the LCIO file with the alignment constants (add .slcio)",
                                _alignmentDBFileName, static_cast< string > ( "eudet-alignmentDB.slcio" ) ); 
    
    std::vector<int> initIgnoreIDVec;
    registerProcessorParameter ("IgnoreIDs",
                                "Ignore clusters from list of sensorIDs",
                                _ignoreIDVec, initIgnoreIDVec);
    
  }
  
  //
  // Method called at the beginning of data processing
  //
  void ClusterShapeToAsciiPrinter::init() {
    
    // Initialize variables
    _nRun = 0 ;
    _nEvt = 0 ;
    _timeCPU = clock()/1000 ;
    
    // Print set parameters
    printProcessorParams();
    
    // Read detector constants from gear file
    _detector.ReadGearConfiguration();  
    
    // Read alignment data base file 
    _detector.ReadAlignmentDB( _alignmentDBFileName );     

    // Open file stream
    _outfile.open(_asciiFileName, ios::out | ios::trunc);
    _outfile << "sensor\t"
             << "label\t"
             << "trk_clu_u[mm]\t"
             << "trk_clu_v[mm]\t"
             << "trk_dudw\t"
             << "trk_dvdw\t"
             << "trk_mom[GeV]\t"
             << "trk_cov_uu[mm2]\t"
             << "trk_cov_vv[mm2]\t"
             << "trk_cov_uv[mm2]"
             << endl;
  }
  
  //
  // Method called for each run
  //
  void ClusterShapeToAsciiPrinter::processRunHeader(LCRunHeader * run)
  {
   
    // Print run number
    streamlog_out(MESSAGE3) << "Processing run: "
                            << (run->getRunNumber())
                            << std::endl << std::endl;
    
    _nRun++ ;

  }
  
  //
  // Method called for each event
  //
  void ClusterShapeToAsciiPrinter::processEvent(LCEvent * evt)
  {
    
    _nEvt ++ ;
    
    TrackInputProvider TrackIO; 
    
    GenericTrackFitter TrackFitter(_detector);
    TrackFitter.SetNumIterations(2); 
    
    LCCollection* inputCollection;
    try {
      inputCollection = evt->getCollection(_inputTrackCollectionName);
    } catch (DataNotAvailableException& e) {
      throw SkipEventException(this);
    }
     
    // Main loop over all tracks
    int nTracks = inputCollection->getNumberOfElements(); 
    for (int itrk = 0; itrk < nTracks; itrk++) {   
      
      // Retrieve track from LCIO 
      Track * inputtrack = dynamic_cast<Track*> (inputCollection->getElementAt(itrk));
      
      // Convert LCIO -> TB track  
      TBTrack track = TrackIO.MakeTBTrack( inputtrack, _detector );  
      
      // Refit track 
      bool trkerr = TrackFitter.Fit(track);
      if ( trkerr ) {
        continue;
      } 
      
      //
      // Loop over all clusters in the track
      for (int ipl= 0; ipl< _detector.GetNSensors(); ++ipl) {   
        
        // Get sensor data 
        //------------------------
        TBTrackElement& TE = track.GetTE(ipl);  
        Det & Sensor = _detector.GetDet(ipl);  
        int sensorID = Sensor.GetDAQID();        
        
        bool ignoreID = false;
        for (auto id :  _ignoreIDVec)  {
          if  (id == sensorID) ignoreID = true; 
        }
         
        // Ignore track elements w/o measurment
        if ( TE.HasHit() && !ignoreID ) { 
          
          // Get local track parameters 
          double trk_tu = TE.GetState().GetPars()[0][0];  // rad
          double trk_tv = TE.GetState().GetPars()[1][0];  // rad
          double trk_u = TE.GetState().GetPars()[2][0];   // mm
          double trk_v = TE.GetState().GetPars()[3][0];   // mm
          double trk_qp = TE.GetState().GetPars()[4][0];  // 1/GeV
          double trk_charge = track.GetCharge();
          double trk_mom = std::abs(trk_charge/trk_qp); 
           
          double sigma2_u = TE.GetState().GetCov()[2][2]; 
          double sigma2_v = TE.GetState().GetCov()[3][3]; 
          double cov_uv   = TE.GetState().GetCov()[2][3];  
           
          PixelCluster Cluster = TE.GetHit().GetCluster();  
          string id = Cluster.getLabel();
          
          trk_u -= Sensor.GetPixelCenterCoordU( Cluster.getVStart(), Cluster.getUStart()); 
          trk_v -= Sensor.GetPixelCenterCoordV( Cluster.getVStart(), Cluster.getUStart()); 
           
          // Dump these numbers into ascii file
          _outfile << ipl << "\t"
                   << id.c_str() << "\t"
                   << trk_u << "\t"
                   << trk_v << "\t"
                   << trk_tu << "\t"
                   << trk_tv << "\t"
                   << trk_mom << "\t"
                   << sigma2_u << "\t"
                   << sigma2_v << "\t"
                   << cov_uv
                   << endl;
           
        }
      }
    }  
    
    return;
  }
  
  
  //
  // Method called after each event to check the data processed
  //
  void ClusterShapeToAsciiPrinter::check( LCEvent * evt ) {}
  
  //
  // Method called after all data processing
  //
  void ClusterShapeToAsciiPrinter::end()
  {
    
    // CPU time end
    _timeCPU = clock()/1000 - _timeCPU;
    
    // Print message
    streamlog_out(MESSAGE3) << std::endl
                            << " "
                            << "Time per event: "
                            << std::setiosflags(std::ios::fixed | std::ios::internal )
                            << std::setprecision(3)
                            << _timeCPU/_nEvt
                            << " ms"
                            << std::endl
                            << std::setprecision(3)
                            << std::endl
                            << " "
                            << "Processor succesfully finished!"
                            << std::endl;
     
    
    // Close file stream
    _outfile << endl;
    _outfile.close();
  }

  //
  // Method printing processor parameters
  //
  void ClusterShapeToAsciiPrinter::printProcessorParams() const 
  {
    
    streamlog_out(MESSAGE3)  << std::endl
                             << " "
                             << "ClusterShapeToAsciiPrinter Development Version, be carefull!!"
                             << " "
                             << std::endl  << std::endl;   
    
  }
  
  
} // Namespace



