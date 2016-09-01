// EventViewer Processor  
// 
// See EventViewer.h for full documentation of processor. 
// 
// Author: Benjamin Schwenker, Göttingen University 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "EventViewer.h"

// Include basic C
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

// Include LCIO classes
#include <lcio.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// Include TBTools 
#include "Utilities.h"


// Used namespaces
using namespace std; 
using namespace lcio ;
using namespace marlin ;

namespace depfet {

//
// Instantiate this object
//
EventViewer aEventViewer ;

//
// Constructor
//
EventViewer::EventViewer() : Processor("EventViewer")
{

// Processor description
   _description = "EventViewer: Depfet event display" ;

//
// Processor parameters

   // Define compulsory parameters
   registerInputCollection(    LCIO::TRACKERDATA,
                               "PixelDataCollection" ,
                               "Name of pixel data collection"  ,
                               _DataCollectionName,
                               std::string("data") ) ;

   registerInputCollection(    LCIO::TRACKERRAWDATA, 
                               "PixelRawDataCollection",
			       "Name of pixel raw data collection",
			       _RawDataCollectionName, 
                               string ("rawdata"));
   
   
   registerProcessorParameter( "RootFileName",
                               "Output root file name",
                               _rootFileName,
                               std::string("TB"));
   
   registerProcessorParameter("SelectDataType", "Available types are: \n"
                             "RAWDATA: NxM matrix of raw data\n"
                             "DATA: NxM matrix of corr. data \n"
                             "ZEROSUPP: Zero suppressed data",
                             _selectDataType, string( "RAWDATA" ) );
   
   registerProcessorParameter("SelectTriggerType", "Available types are: \n"
                             "EXTERNAL: Use list of triggered events from trigger file\n"
                             "ALL: Trigger all events\n",
                             _selectTrigger, string( "ALL" ) );
   
   registerProcessorParameter( "MaxTriggers",
                               "Max number of triggers",
                               _maxTriggers,
                               static_cast<int>(1000));
   
   registerProcessorParameter( "UseStatusMap",
                               "Display only GOOD Pixels, not available for RAWDATA",
                               _useStatusMap,
                               static_cast<bool>(false));
   
   registerProcessorParameter( "DisplaySensorID",
                               "-1: for all sensors or use sensorID to view specific sensor",
                               _displaySensorID,
                               static_cast<int>(-1));

   registerProcessorParameter( "TiggerFileName",
                               "Input file with triggered event numbers",
                               _triggerFileName,
                               std::string("dummy.txt"));
  
}

//
// Method called at the beginning of data processing
//
void EventViewer::init() {

// Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   _triggerCounter = 0;
   
// Print set parameters
   printProcessorParams();

// Read detector constants from gear file
   _detector.ReadGearConfiguration();    
   
// transform the selectDataType string to small letters
   transform(_selectDataType.begin(), _selectDataType.end(), _selectDataType.begin(), ::tolower);
   transform(_selectTrigger.begin(), _selectTrigger.end(), _selectTrigger.begin(), ::tolower);

// CPU time start
   _timeCPU = clock()/1000;

// ROOT_OUTPUT
   string hname = _rootFileName + "_DepfetEvent.root"; 
   _rootFile = new TFile(hname.c_str(),"recreate");

// Check if trigger file needed
   _triggers.clear();  
   if ( _selectTrigger == "external" ) {
          
     string line;
     ifstream triggerFile( _triggerFileName.c_str() );
     if ( triggerFile.is_open() ) { 
       
       streamlog_out(MESSAGE3) << "Reading trigger file ... "<< endl; 
       while (! triggerFile.eof() )
       {
         getline ( triggerFile ,line);
         streamlog_out(MESSAGE2) << "reading line: "<< line; 

         int i = -1;
         
         if( sscanf( line.c_str() , "%d", &i)  == 1 && i>-1 )
         {
           streamlog_out(MESSAGE2) << "  converted eventID "<< i << endl;  
           _triggers.push_back(i);
         } else {
           streamlog_out(MESSAGE2) << "  bad eventID, skip " << endl;  
         }

       }
       triggerFile.close();

     } else {
       streamlog_out(ERROR) << "Unable to open trigger file. " << endl; 
       exit(1); 
     }
     
     // sort triggers to improve search speed
     sort (_triggers.begin(), _triggers.end());

   } // endif 
   
}

//
// Method called for each run
//
void EventViewer::processRunHeader(LCRunHeader * run)
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
void EventViewer::processEvent(LCEvent * evt)
{
   
   
 
   // Print event number
   if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                   << (evt->getEventNumber())
                                                                   << std::endl << std::endl;
   
   if ( _triggerCounter > _maxTriggers ) {
     streamlog_out( MESSAGE3 ) << "Enough triggers processed!" << endl;
     throw StopProcessingException(this);
   }
           
   //
   // Open collections
   try {
      
     // check if event is triggered
     bool trigger = false; 
     if ( _selectTrigger == "external" ) {        
         trigger = checkExternalTrigger(evt);    
     } else  {
         trigger = true; 
     }
        
     if ( trigger ) {
         streamlog_out(MESSAGE3) << "Dumping event "<< evt->getEventNumber() <<" as number "<<_triggerCounter<<endl;
         if  ( _selectDataType == "data" ) dumpDataEvent( evt );
         else if ( _selectDataType  == "rawdata" )  dumpRawDataEvent( evt );
         else if ( _selectDataType  == "zerosupp" )  dumpZeroSuppEvent( evt );
         else {
           exit(-1); 
         }
         _triggerCounter++;          
     }
       
   } catch(DataNotAvailableException &e){
      streamlog_out(ERROR4) << "Continuing with next event" << std::endl;
   }  
   _nEvt ++ ;
}

//
// Looking for eventID in list of selected events
// 
bool EventViewer::checkExternalTrigger( LCEvent * evt )
{
 
  int eventID = evt->getEventNumber();
  bool trigger = binary_search (_triggers.begin(), _triggers.end(), eventID); 
  return trigger; 

}

//
// Method dumping event data frames to 2D histos
//
void EventViewer::dumpDataEvent( LCEvent * evt )
{
  
   
  int eventID = evt->getEventNumber();  
  
  //
  // Open collections
  try {
    
    LCCollectionVec * frames = dynamic_cast < LCCollectionVec * > (evt->getCollection(_DataCollectionName)); 
      
    if ( !_useStatusMap ) {    
      
      // Loop over all sensors 
      for (unsigned int iSensor = 0; iSensor < frames->size(); iSensor++) {
        
        //
        // open data frame
        TrackerDataImpl * matrix = dynamic_cast<TrackerDataImpl* > (frames->getElementAt(iSensor));
        CellIDDecoder<TrackerDataImpl> idMatrixDecoder(frames);
        int currentSensorID = idMatrixDecoder(matrix)["sensorID"]; 
        int noOfXPixels =  idMatrixDecoder(matrix)["xMax"]+1;    
        int noOfYPixels =  idMatrixDecoder(matrix)["yMax"]+1;   

        //
        // skip this sensor 
        if ( _displaySensorID != currentSensorID  && _displaySensorID != -1 ) continue;     
 
        streamlog_out(MESSAGE4) << "Dump sensor data wo status "<<currentSensorID<<endl;  

        // 
        // prepare histo 
        std::string histoName = _rootFileName+"_evt_"+to_string(eventID)+"_mod_"+to_string(currentSensorID); 
        std::string  histoTitle = "Evt:"+to_string(eventID)+" Mod:"+to_string(currentSensorID);
      	   
        double xmin = 0 - 0.5;
        double xmax = noOfXPixels -1 + 0.5;
        int xnbins = noOfXPixels;
        double ymin = 0 - 0.5;
        double ymax = noOfYPixels -1 + 0.5;
        int ynbins = noOfYPixels; 
        TH2D * eventMap = new TH2D(histoName.c_str(),histoTitle.c_str(),xnbins,xmin,xmax,ynbins,ymin,ymax);
        eventMap->SetXTitle("X Axis"); 
        eventMap->SetYTitle("Y Axis");

        //
        // fill 2D histo
        FloatVec charges  = matrix->getChargeValues();    
        int iPixel = 0; 
        for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
          for ( int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
        
            double chargeValue = charges[iPixel]; 
            eventMap->Fill(xPixel,yPixel,chargeValue); 
            
            iPixel++;  
      
          }
        }
     
        // write map to root file 
        _rootFile->cd("/"); 
        eventMap->Write();
    
      } // end sensor loop
        
    } else {
      
      // Load status collection   
      LCCollectionVec * statusCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection("status"));  

      // Loop over all sensors 
      for (unsigned int iSensor = 0; iSensor < frames->size(); iSensor++) {
       
        //
        // open data frame
        TrackerDataImpl * matrix = dynamic_cast<TrackerDataImpl* > (frames->getElementAt(iSensor));
        CellIDDecoder<TrackerDataImpl> idMatrixDecoder(frames);
        int currentSensorID = idMatrixDecoder(matrix)["sensorID"]; 
        int noOfXPixels =  idMatrixDecoder(matrix)["xMax"]+1;    
        int noOfYPixels =  idMatrixDecoder(matrix)["yMax"]+1;   

        // open status map 
        TrackerRawDataImpl * status = dynamic_cast<TrackerRawDataImpl* > (statusCollection->getElementAt(iSensor));
        
        //
        // skip this sensor 
        if ( _displaySensorID != currentSensorID  && _displaySensorID != -1 ) continue;     
 
        streamlog_out(MESSAGE1) << "Dump sensor "<<currentSensorID<<endl;  

        // 
        // prepare histo 
        std::string histoName = _rootFileName+"_evt_"+to_string(eventID)+"_mod_"+to_string(currentSensorID); 
        std::string  histoTitle = "Evt:"+to_string(eventID)+" Mod:"+to_string(currentSensorID);
      	   
        double xmin = 0 - 0.5;
        double xmax = noOfXPixels -1 + 0.5;
        int xnbins = noOfXPixels;
        double ymin = 0 - 0.5;
        double ymax = noOfYPixels -1 + 0.5;
        int ynbins = noOfYPixels; 
        TH2D * eventMap = new TH2D(histoName.c_str(),histoTitle.c_str(),xnbins,xmin,xmax,ynbins,ymin,ymax);
        eventMap->SetXTitle("X Axis"); 
        eventMap->SetYTitle("Y Axis");

        //
        // fill 2D histo
        FloatVec charges  = matrix->getChargeValues();
        ShortVec pixelStatus = status->getADCValues();     
        int iPixel = 0; 
        for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
          for ( int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
        
            double chargeValue = charges[iPixel];
            short statusValue = pixelStatus[iPixel]; 
         
            if ( statusValue != 0 ) {
              eventMap->Fill(xPixel,yPixel,0.0); 
            } else {
              eventMap->Fill(xPixel,yPixel,chargeValue); 
            }

            iPixel++;  
      
          }
        }
        
        // write map to root file 
        _rootFile->cd("/"); 
        eventMap->Write();
    
      } // end sensor loop
   
    } // end _useStatus

  } catch(DataNotAvailableException &e){
     streamlog_out(ERROR4) << "No data and/or status collection found" << std::endl;
  }  

  
}

//
// Method dumping event rawdata frames to 2D histos
//
void EventViewer::dumpRawDataEvent( LCEvent * evt )
{

  int eventID = evt->getEventNumber();  
  
  //
  // Open collections
  try {
    
    LCCollectionVec * frames = dynamic_cast < LCCollectionVec * > (evt->getCollection(_RawDataCollectionName)); 
            
    // Loop over all sensors 
    for (unsigned int iSensor = 0; iSensor < frames->size(); iSensor++) {

      //
      // open data frame
      TrackerRawDataImpl * rawmatrix = dynamic_cast<TrackerRawDataImpl* > (frames->getElementAt(iSensor));
      CellIDDecoder<TrackerRawDataImpl> idRawMatrixDecoder( frames );
      int currentSensorID = idRawMatrixDecoder(rawmatrix)["sensorID"]; 
      int noOfXPixels =  idRawMatrixDecoder(rawmatrix)["xMax"]+1;    
      int noOfYPixels =  idRawMatrixDecoder(rawmatrix)["yMax"]+1;   

      //
      // skip this sensor 
      if ( _displaySensorID != currentSensorID  && _displaySensorID != -1 ) continue;     
 
      streamlog_out(MESSAGE1) << "Dump sensor "<<currentSensorID<<endl;  

      // 
      // prepare histo 
      std::string histoName = _rootFileName+"_rawEvt_"+to_string(eventID)+"_mod_"+to_string(currentSensorID); 
      std::string  histoTitle = "Raw Evt:"+to_string(eventID)+" Mod:"+to_string(currentSensorID);
      	   
      double xmin = 0 - 0.5;
      double xmax = noOfXPixels -1 + 0.5;
      int xnbins = noOfXPixels;
      double ymin = 0 - 0.5;
      double ymax = noOfYPixels -1 + 0.5;
      int ynbins = noOfYPixels; 
      TH2D * eventMap = new TH2D(histoName.c_str(),histoTitle.c_str(),xnbins,xmin,xmax,ynbins,ymin,ymax);
      eventMap->SetXTitle("X Axis"); 
      eventMap->SetYTitle("Y Axis");

      //
      // fill 2D histo    
      ShortVec charges = rawmatrix->getADCValues();
      int iPixel = 0; 
      for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
        for ( int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
        
          double chargeValue = charges[iPixel];   
          eventMap->Fill(xPixel,yPixel,chargeValue);   
          iPixel++;  
      
        }
      }
     
      // write map to root file 
      _rootFile->cd("/"); 
      eventMap->Write();

    } // end sensor loop

  } catch(DataNotAvailableException &e){
     streamlog_out(ERROR4) << "No raw data collection found" << std::endl;
  }   
  
}

//
// Method dumping event data from Mimosa26 
//
void EventViewer::dumpZeroSuppEvent( LCEvent * evt )
{

  int eventID = evt->getEventNumber();   

  

  //
  // Open collections
  try {
    
    LCCollectionVec * frames = dynamic_cast < LCCollectionVec * > (evt->getCollection(_DataCollectionName));           
    
         
    
    // Loop over all sensors 
    for (unsigned int iSensor = 0; iSensor < frames->size(); iSensor++) {
      
     
      CellIDDecoder<TrackerDataImpl> idMatrixDecoder( frames );
      
      //
      // open data frame
      TrackerDataImpl * matrix = dynamic_cast<TrackerDataImpl* > (frames->getElementAt(iSensor));
      int currentSensorID = idMatrixDecoder(matrix)["sensorID"]; 
      int ipl = _detector.GetPlaneNumber(currentSensorID);   
           
      streamlog_out(MESSAGE4) << "Dump sensor iDetector " << iSensor << " with sensorID " << currentSensorID << endl;  
      
      int noOfXPixels =  _detector.GetDet(ipl).GetNColumns();   
      int noOfYPixels =  _detector.GetDet(ipl).GetNRows();  
      
      //
      // skip this sensor 
      if ( _displaySensorID != currentSensorID  && _displaySensorID != -1 ) continue;     
      
      double xmin = 0 - 0.5;
      double xmax = noOfXPixels -1 + 0.5;
      int xnbins = noOfXPixels;
      double ymin = 0 - 0.5;
      double ymax = noOfYPixels -1 + 0.5;
      int ynbins = noOfYPixels; 
      	   
      TH2D * eventMap = new TH2D(Form("evt%d_mod%d",eventID,currentSensorID),Form("Data Event %d",eventID),xnbins,xmin,xmax,ynbins,ymin,ymax);
      eventMap->SetXTitle("columns"); 
      eventMap->SetYTitle("rows");

      
      FloatVec sparsePixels = matrix->getChargeValues();
      int nPixels = sparsePixels.size()/3; 
      
      //
      // fill 2D histo   
      for ( int index=0; index<nPixels;  index++) { 

        int xPixel = static_cast<int> (sparsePixels[index * 3]);
        int yPixel = static_cast<int> (sparsePixels[index * 3 + 1]);
        float chargeValue =  sparsePixels[index * 3 + 2];
       
         
        eventMap->Fill(xPixel,yPixel,chargeValue); 
        
      }
       
      // write map to root file 
      _rootFile->cd("/"); 
      eventMap->Write();
      

    } // end sensor loop

  } catch(DataNotAvailableException &e){
     streamlog_out(ERROR4) << "No raw data collection found" << std::endl;
  }   
  
}


//
// Method called after each event to check the data processed
//
void EventViewer::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void EventViewer::end()
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

    
   // Close file
   _rootFile->Write();
   _rootFile->Close();

}


//
// Method printing processor parameters
//
void EventViewer::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "EventViewer Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

} // Namespace
