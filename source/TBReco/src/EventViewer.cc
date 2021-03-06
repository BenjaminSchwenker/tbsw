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
#include "TBDetector.h"


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
                               std::string("EventViewer.root"));
   
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
   
   
   
   // transform the selectDataType string to small letters
   transform(_selectDataType.begin(), _selectDataType.end(), _selectDataType.begin(), ::tolower);
   transform(_selectTrigger.begin(), _selectTrigger.end(), _selectTrigger.begin(), ::tolower);
   
   // CPU time start
   _timeCPU = clock()/1000;
   
   // ROOT_OUTPUT
   _rootFile = new TFile(_rootFileName.c_str(),"recreate");
   
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
       
       // sort triggers to improve search speed
       sort (_triggers.begin(), _triggers.end());
       
     } else {
       streamlog_out(ERROR) << "Unable to open trigger file. " << endl;  
     }
     
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
   
             
   //
   // Open collections
   try {
      
     // check if event is triggered
     bool trigger = false; 
     if ( _selectTrigger == "external" ) {        
         trigger = checkExternalTrigger(evt);    
     } else  {
         if ( _triggerCounter < _maxTriggers ) trigger = true; 
     }
        
     if ( trigger ) {
         streamlog_out(MESSAGE3) << "Dumping event "<< evt->getEventNumber() <<" as number "<<_triggerCounter<<endl;
         if  ( _selectDataType == "data" ) dumpDataEvent( evt );
         else if ( _selectDataType  == "rawdata" )  dumpRawDataEvent( evt );
         else if ( _selectDataType  == "zerosupp" )  dumpZeroSuppEvent( evt );
         else {
           streamlog_out(ERROR) << "Invalid data type selected. Should be 'data' or 'rawdata' or 'zerosupp'." << std::endl;
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
        
    // Loop over all sensors 
    for (unsigned int iSensor = 0; iSensor < frames->size(); iSensor++) {
        
        //
        // open data frame
        TrackerDataImpl * matrix = dynamic_cast<TrackerDataImpl* > (frames->getElementAt(iSensor));
        CellIDDecoder<TrackerDataImpl> idMatrixDecoder(frames);
        int currentSensorID = idMatrixDecoder(matrix)["sensorID"];

        //
        // skip this sensor 
        if ( _displaySensorID != currentSensorID  && _displaySensorID != -1 ) continue;     
 
        streamlog_out(MESSAGE3) << "Dump sensor data wo status "<<currentSensorID<<endl;  
 
        int ipl = TBDetector::GetInstance().GetPlaneNumber(currentSensorID);   
        int noOfXPixels =  TBDetector::Get(ipl).GetMaxUCell()-TBDetector::Get(ipl).GetMinUCell()+1;   
        int noOfYPixels =  TBDetector::Get(ipl).GetMaxVCell()-TBDetector::Get(ipl).GetMinVCell()+1; 
       
        // 
        // create histo in local root file
        _rootFile->cd("");  
        std::string histoName = _rootFileName+"_evt_"+to_string(eventID)+"_mod_"+to_string(currentSensorID); 
        std::string  histoTitle = "Evt:"+to_string(eventID)+" Mod:"+to_string(currentSensorID);
      	   
        double xmin = TBDetector::Get(ipl).GetMinUCell();
        double xmax = noOfXPixels + xmin;
        int xnbins = noOfXPixels;
        double ymin = TBDetector::Get(ipl).GetMinVCell();
        double ymax = noOfYPixels + ymin;
        int ynbins = noOfYPixels; 
        TH2D * eventMap = new TH2D(histoName.c_str(),histoTitle.c_str(),xnbins,xmin-0.5,xmax-0.5,ynbins,ymin-0.5,ymax-0.5);
        eventMap->SetXTitle("X Axis"); 
        eventMap->SetYTitle("Y Axis");

        //
        // fill 2D histo
        FloatVec charges  = matrix->getChargeValues();    
        int iPixel = 0; 
        for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
          for ( int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
        
            double chargeValue = charges[iPixel]; 
            eventMap->Fill(xPixel+xmin,yPixel+ymin,chargeValue); 
            
            iPixel++;  
      
          }
        }
     
        // write map to local root file 
        eventMap->Write();
    
    } // end sensor loop
        
  } catch(DataNotAvailableException &e){
     streamlog_out(ERROR4) << "No data collection found" << std::endl;
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
      int ipl = TBDetector::GetInstance().GetPlaneNumber(currentSensorID);   

      //
      // skip this sensor 
      if ( _displaySensorID != currentSensorID  && _displaySensorID != -1 ) continue;     
      
      streamlog_out(MESSAGE3) << "Dump sensor plane " << ipl << " with sensorID " << currentSensorID << endl;  
 
      int noOfXPixels =  TBDetector::Get(ipl).GetMaxUCell()-TBDetector::Get(ipl).GetMinUCell()+1;   
      int noOfYPixels =  TBDetector::Get(ipl).GetMaxVCell()-TBDetector::Get(ipl).GetMinVCell()+1;  
     
      // 
      // prepare histo 
      _rootFile->cd(""); 
      std::string histoName = _rootFileName+"_rawEvt_"+to_string(eventID)+"_mod_"+to_string(currentSensorID); 
      std::string  histoTitle = "Raw Evt:"+to_string(eventID)+" Mod:"+to_string(currentSensorID);
      	   
      double xmin = TBDetector::Get(ipl).GetMinUCell();
      double xmax = noOfXPixels + xmin;
      int xnbins = noOfXPixels;
      double ymin = TBDetector::Get(ipl).GetMinVCell();
      double ymax = noOfYPixels + ymin;
      int ynbins = noOfYPixels; 
      TH2D * eventMap = new TH2D(histoName.c_str(),histoTitle.c_str(),xnbins,xmin-0.5,xmax-0.5,ynbins,ymin-0.5,ymax-0.5);
      eventMap->SetXTitle("X Axis"); 
      eventMap->SetYTitle("Y Axis");

      // 
      // fill 2D histo    
      ShortVec charges = rawmatrix->getADCValues();
      int iPixel = 0; 
      for (int yPixel = 0; yPixel < noOfYPixels; yPixel++) {
        for ( int xPixel = 0; xPixel < noOfXPixels; xPixel++) {
        
          double chargeValue = charges[iPixel];   
          eventMap->Fill(xPixel+xmin,yPixel+ymin,chargeValue);   
          iPixel++;  
      
        }
      }
     
      // write map to root file 
      eventMap->Write();

    } // end sensor loop

  } catch(DataNotAvailableException &e){
     streamlog_out(ERROR4) << "No raw data collection found" << std::endl;
  }   
  
}

//
// Method dumping zero suppressed event data 
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
      
      //
      // skip this sensor 
      if ( _displaySensorID != currentSensorID  && _displaySensorID != -1 ) continue;     

      streamlog_out(MESSAGE3) << "Dump sensor iDetector " << iSensor << " with sensorID " << currentSensorID << endl;  

      int ipl = TBDetector::GetInstance().GetPlaneNumber(currentSensorID);   
      int noOfXPixels =  TBDetector::Get(ipl).GetMaxUCell()-TBDetector::Get(ipl).GetMinUCell()+1;   
      int noOfYPixels =  TBDetector::Get(ipl).GetMaxVCell()-TBDetector::Get(ipl).GetMinVCell()+1;  
     
      // 
      // prepare histo 
      _rootFile->cd(""); 
      double xmin = TBDetector::Get(ipl).GetMinUCell();
      double xmax = noOfXPixels + xmin;
      int xnbins = noOfXPixels;
      double ymin = TBDetector::Get(ipl).GetMinVCell();
      double ymax = noOfYPixels + ymin;
      int ynbins = noOfYPixels; 
   
      TH2D * eventMap = new TH2D(Form("evt%d_mod%d",eventID,currentSensorID),Form("Data Event %d",eventID),xnbins,xmin-0.5,xmax-0.5,ynbins,ymin-0.5,ymax-0.5);
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
      eventMap->Write();
      

    } // end sensor loop

  } catch(DataNotAvailableException &e){
     streamlog_out(ERROR4) << "No zerosupp data collection found" << std::endl;
  }   
  
}


//
// Method called after each event to check the data processed
//
void EventViewer::check( LCEvent * )
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
