// StripClusterizer implementation file
// 
// Author: Benjamin Schwenker, University of Göttingen 
// <mailto:benjamin.schwenker@phys.uni-goettingen.de>

#include "StripClusterizer.h"

// Include DEPFETTrackTools 
#include "DEPFET.h" 

#include <iomanip>

// Used namespaces
using namespace std; 
using namespace lcio;
using namespace marlin;

namespace depfet {

//
// Instantiate this object
//
StripClusterizer aStripClusterizer ;

//
// Constructor
//
StripClusterizer::StripClusterizer() : Processor("StripClusterizer")
{

// Processor description
   _description = "StripClusterizer: Looking for clusters in strip digits" ;

//   
// First of all, we need to register the input/output collections
   
   registerInputCollection (LCIO::TRACKERDATA, "SparseDataCollectionName",
                            "Name of input sparsified strip digit collection",
                            _sparseDataCollectionName, string("zsdata"));
   
   registerInputCollection (LCIO::TRACKERRAWDATA, "StatusCollectionName",
                           "Input status collection name",
                           _statusCollectionName, string("status")); 
    
   registerOutputCollection (LCIO::TRACKERPULSE, "ClusterCollectionName",
                            "Name of the output cluster collection",
                            _clusterCollectionName, string("zscluster"));
    
   registerProcessorParameter( "SparseZSCut","Threshold for zero suppression",
                               _sparseZSCut, static_cast<float > (0));

   registerProcessorParameter( "SparseSeedCut","Threshold for seed digits",
                               _sparseSeedCut, static_cast<float > (5));
   
   registerProcessorParameter( "SparseClusterCut","Threshold for cluster signal",
                               _sparseClusterCut, static_cast<float > (7) ); 
  
   registerProcessorParameter( "AcceptGaps","Accepts clusters with N missing strips",
                               m_acceptGaps , static_cast<int > (1) ); 

   registerProcessorParameter( "nSamples","Number of strip signal samples per hit",
                               m_samples , static_cast<int > (1) ); 
    
}

//
// Method called at the beginning of data processing
//
void StripClusterizer::init() {
   
   // Initialize variables
   _nRun = 0 ;
   _nEvt = 0 ;
   
   // Read detector constants from gear file
   _detector.ReadGearConfiguration();    
   
   _dummyCollectionName = "original_data_"+_clusterCollectionName;
   
   // The status data is not yet initialized
   _isStatusReady = false;
              
   // Print set parameters
   printProcessorParams();
   
   // CPU time start
   _timeCPU = clock()/1000;
   
}

//
// Method called for each run
//
void StripClusterizer::processRunHeader(LCRunHeader * run)
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
void StripClusterizer::processEvent(LCEvent * evt)
{
   
   // Print event number
   if ((evt->getEventNumber())%100 == 0) streamlog_out(MESSAGE3) << "Events processed: "
                                                                   << (evt->getEventNumber())
                                                                   << std::endl << std::endl;
   
   // More detailed event numbering for testing
   streamlog_out(MESSAGE2) << std::endl << "Starting with Event Number " << evt->getEventNumber()  << std::endl;  
   
   // First of all we need to be sure that status data is initialized
   if ( !_isStatusReady ) {
     initializeStatus( evt ) ;
   }
   
   //
   // Open collections
   try {
     
     // Output collection containing cluster data 
     LCCollectionVec * clusterCollection = new LCCollectionVec(LCIO::TRACKERPULSE) ;
     
     // Build clusters, i.e. adjecant groups of sparsified digits passing certain signal cuts  
     clusterize( evt , clusterCollection  ); 
     
     // Add clusterCollection to event
     evt->addCollection( clusterCollection, _clusterCollectionName );
             	    
   } catch(DataNotAvailableException &e){}  
   _nEvt ++ ;
}


//
// Method called after each event to check the data processed
//
void StripClusterizer::check( LCEvent * evt )
{
}

//
// Method called after all data processing
//
void StripClusterizer::end()
{
   
   // CPU time end
   _timeCPU = clock()/1000 - _timeCPU;
   
   // Print message
   streamlog_out(MESSAGE3) << std::endl
                           << " "
                           << "Time per event: "
                           << std::setprecision(3)
                           << _timeCPU/_nEvt
                           << " ms"
                           << std::endl
                           << std::setprecision(3)
                           << std::endl
                           << " "
                           << "Processor succesfully finished!"
                           << std::endl;
   
   
   
}


//
// Method printing processor parameters
//
void StripClusterizer::printProcessorParams() const
{

   streamlog_out(MESSAGE3)  << std::endl
                            << " "
                            << "StripClusterizer Development Version, be carefull!!"
                            << " "
                            << std::endl  << std::endl;   


}

// Called by the processEvent() 

void StripClusterizer::clusterize( LCEvent * evt , LCCollectionVec * clusterCollection  ) 
{  
    
  // Open zero suppressed digits 
  LCCollectionVec * DigitCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_sparseDataCollectionName)); 
  // Helper class for decoding strip data 
  CellIDDecoder<TrackerDataImpl> StripID( DigitCollection );  
  
  //Open status data
  LCCollectionVec * StatusCollection = dynamic_cast < LCCollectionVec * > (evt->getCollection(_statusCollectionName));
  
  // The original data collection contains the digits for all clusters
  LCCollectionVec * originalDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);
  
  // Clustering loop over strip detectors 
  for (unsigned int iDet = 0; iDet < DigitCollection->size(); iDet++) { 
     
    // Get strip digits from next detector   
    TrackerDataImpl * Digits = dynamic_cast<TrackerDataImpl* > ( DigitCollection->getElementAt(iDet) );
     
    // Get status matrix 
    TrackerRawDataImpl * Status = dynamic_cast<TrackerRawDataImpl*> (StatusCollection->getElementAt( iDet ));  
      
    // Get DAQ ID 
    int sensorID = StripID( Digits ) ["sensorID"];

    // Read geometry info for sensor 
    int ipl = _detector.GetPlaneNumber(sensorID);      
    Det& Sensor = _detector.GetDet(ipl);
  
    int noOfXPixels = Sensor.GetNColumns(); 
    int noOfYPixels = Sensor.GetNRows();

    // Get max channel numbers 
    int maxCol = Sensor.GetNColumns() - 1; 
    int maxRow = Sensor.GetNRows() - 1;
    
    
    FloatVec rawDigits = Digits->getChargeValues();
    int nDigits = rawDigits.size()/(3*m_samples); 
        
    // Group all digits into cluster candidates, groups of geometrically adjecant 
    // or neighboring digits. 
    
    ClusterCandVec Clusters;
    Clusters.reserve(nDigits);
      
    // Loop over digits
    for (int i = 0; i < nDigits; i++) 
    {   
      
      // Read first sample of hit measurement
      int isU = static_cast<int> (rawDigits[i * 3 * m_samples]);
      int cell = static_cast<int> (rawDigits[i * 3 * m_samples + 1]);
      float signal = rawDigits[i * 3 * m_samples + 2]; 
      
      // Loop over remaining six samples and find max. signal 
      for (int is = 1; is < m_samples; is++) {
        float next_signal = rawDigits[i * 3 * m_samples + is * 3 + 2]; 
        if ( next_signal > signal ) signal = next_signal; 
      }

      streamlog_out(MESSAGE1) << "Digits " << i << " on sensor " << sensorID  
                              << std::endl;  
      streamlog_out(MESSAGE1) << "   cell: " << cell << ", isU: " << isU
                              << ", charge: " << signal
                              << std::endl;

      // If signal below threshold, skip it in clusterization
      if ( signal < _sparseZSCut ) {
        streamlog_out(MESSAGE2) << "  Signal below ZS cut. Skipping it." << std::endl; 
        continue;
      }
        
      /*
      // If a pixel is out of range, skip it in clusterization
      if ( col < 0 || col > maxCol || row < 0 || row > maxRow ) {
        streamlog_out(MESSAGE2) << "  Invalid pixel address found. Skipping it." << std::endl; 
        continue;
      }
      
      // If a pixel is bad, skip it in clusterization
      if ( status->getADCValues() [ matrixDecoder.getIndexFromXY(col,row) ]   != 0 ) {
        streamlog_out(MESSAGE2) << "  Bad pixel found. Skipping it." << std::endl; 
        continue;
      }
      */
      
      // Loop on all existing clusters, until you find that current digit 
      // is a neighbour of a cluster.
      bool found= false;
         	
      ClusterCandVec::iterator currentCluster = Clusters.begin();
      ClusterCandVec::iterator lastCluster  = Clusters.end();    
              
      while( !found && currentCluster!= lastCluster)
      {
        
        if ( areNeighbours( *currentCluster, cell, isU, m_acceptGaps ) )
        {
           
          // If digit is a duplicate of one in the cluster, do not add it.   
          if(!isDuplicated( *currentCluster, cell , isU)){

           
            // Add digit to this cluster 
            (*currentCluster).push_back(isU);
            (*currentCluster).push_back(cell);
            (*currentCluster).push_back(signal);
            
            // See if cell is a neighbour to any other groups, if yes perform merging 
            checkForMerge(cell, isU, currentCluster, lastCluster);
              
          } else {
            streamlog_out(MESSAGE2) << "  A strip duplicate found. Skipping it." << std::endl; 
          }
          found = true; 
        }
        ++currentCluster;
      }
      
      // If digit is isolated, seed a new candidate cluster. 
      if(!found)
      {
        FloatVec newClu;
        newClu.push_back(isU);
        newClu.push_back(cell);
        newClu.push_back(signal);
        
        Clusters.push_back(newClu);
        
      }
       
    }  
     
    // Count cluster candidates 
    int candNumber = 0;
    
    // Counts clusters
    int clusterID = 0; 
    
    for( ClusterCandVec::iterator cand = Clusters.begin() ;
     	cand!= Clusters.end() ; ++cand) 
    {
      
      candNumber++; 
      
      // If cluster is empty, i.e. it has been merged with another, 
      // do not attempt to make cluster.
      if ((*cand).size()> 0)
      {
          
        // We have found a new cluster candidate. We loop over all digits
        // to calculate S/N ratios. 
        
        int cluQuality = 0;
        float clusterSignal = 0; 
        float seedSignal = 0; 
        
        
        // Prepare a TrackerData to store original data of cluster
        TrackerDataImpl* sparseCluster = new TrackerDataImpl ; 
        
        int nDigits = (*cand).size()/3; 
        
        for ( int i=0; i < nDigits; i++)
        {
  
          int isU = static_cast<int> ( (*cand)[i * 3]);
          int cell = static_cast<int> ( (*cand)[i * 3 + 1]);
          float signal = ( (*cand)[i * 3 + 2]);
                            
          clusterSignal += signal;
         
          if (signal > seedSignal ) {
            seedSignal = signal; 
          }
             
          // Store data in lcio format
          sparseCluster->chargeValues().push_back( isU );
          sparseCluster->chargeValues().push_back( cell );
          sparseCluster->chargeValues().push_back( signal );   
          
        }
            
       
          
        // Verify if the cluster candidates is a good cluster
        if ( ( seedSignal >= _sparseSeedCut ) && 
             ( clusterSignal >= _sparseClusterCut ) ) {
        
                
 	  // New accpeted cluster
          clusterID++; 
          
          // Ok good cluster ... save it 
          CellIDEncoder<TrackerDataImpl> originalDataEncoder( DEPFET::ZSCLUSTERDEFAULTENCODING, originalDataCollection );
          CellIDEncoder<TrackerPulseImpl> clusterEncoder(DEPFET::ZSCLUSTERDEFAULTENCODING, clusterCollection ); 
           
          // Set the ID for this zsCluster
          originalDataEncoder["sensorID"] = sensorID;
          originalDataEncoder["clusterID"] = 0;
          originalDataEncoder["sparsePixelType"] = static_cast<int> (kSimpleSparsePixel);
          originalDataEncoder["quality"] = cluQuality;
          originalDataEncoder.setCellID( sparseCluster );
          originalDataCollection->push_back( sparseCluster );
                             
          TrackerPulseImpl* zsPulse = new TrackerPulseImpl;
          clusterEncoder["sensorID"]  = sensorID;
          clusterEncoder["clusterID"] = 0;
          clusterEncoder["sparsePixelType"] = static_cast<int> (kSimpleSparsePixel);
          clusterEncoder["quality"] = cluQuality;
          clusterEncoder.setCellID( zsPulse );
          
          zsPulse->setCharge( clusterSignal );
          zsPulse->setQuality( cluQuality );
          zsPulse->setTrackerData( sparseCluster );
          clusterCollection->push_back( zsPulse );          
          
           
          
        } else {
          // Cluster candidate thrown away, so clean up memory    
          sparseCluster->chargeValues().clear(); 
          delete sparseCluster;
        }
        
      } // Non empty cluster candidate 
    } // Cluster candidate loop 
    
    // 
    // Free used memory 
     
    for( ClusterCandVec::iterator cand = Clusters.begin() ;
      	cand!= Clusters.end() ; ++cand) 
    {
      (*cand).clear(); 
    }
    Clusters.clear(); 
    
  } // End detector loop 
  
  // Add original data collection to event
  evt->addCollection( originalDataCollection, _dummyCollectionName );
  
}



// Checks if any other cluster candidate neighbours cell. 
// If so, merge with base.  
 
void StripClusterizer::checkForMerge( int cell, int isU,
 ClusterCandVec::iterator baseGroup,
 ClusterCandVec::iterator lastGroup) 
{
  
  // First non-tested group
  ClusterCandVec::iterator nextGroup(baseGroup+1);
   
  for (; nextGroup!= lastGroup; ++nextGroup)
  {              
    if (areNeighbours( *nextGroup, cell, isU, m_acceptGaps ))
    {
      
      int nDigits = (*nextGroup).size()/3; 
       
      for ( int i=0; i < nDigits; i++)      
      {
        float isU1 = (*nextGroup)[i * 3];
        float cell1 = (*nextGroup)[i * 3 + 1];
        float signal1 =  (*nextGroup)[i * 3 + 2];     
        // Copy pixel to base group 
        (*baseGroup).push_back(isU1);
        (*baseGroup).push_back(cell1);
        (*baseGroup).push_back(signal1);
      }
      (*nextGroup).clear();
       
    }
  }
}




bool StripClusterizer::areNeighbours( FloatVec &group, int cell, int isU, int m_acceptGaps ) 
{   
    
  bool match=false;
  int nDigits = group.size()/3; 
  
  for ( int i=0; i < nDigits; i++)
  {
           
    int isU1 = static_cast<int> (group[i * 3]);
    int cell1 = static_cast<int> (group[i * 3 + 1]);
    
    int delta = abs(cell-cell1);
    
    if ( isU1 == isU && delta <= m_acceptGaps ) {
      match = true;
      break;  
    } 
                
  }
    
  return match;
}


                                        
bool StripClusterizer::isDuplicated( FloatVec &group, int cell, int isU) 
{ 
  bool duplicate = false;
  int nDigits = group.size()/3; 
   
  for ( int i=0; i < nDigits; i++)
  {
    int isU1 = static_cast<int> (group[i * 3]);
    int cell1 = static_cast<int> (group[i * 3 + 1]);  
     
    // Duplicate?? 
    if(isU1 == isU && cell1 == cell){
      duplicate = true;
      break; 
    }
  }  
  return duplicate;
}



 
void  StripClusterizer::initializeStatus( LCEvent * event )  {
  
  streamlog_out( MESSAGE3 ) << "Initializing status" << endl;
  
  try {
         
    // Open data  
    LCCollectionVec * DigitCollection = dynamic_cast < LCCollectionVec * > (event->getCollection(_sparseDataCollectionName)); 
    CellIDDecoder<TrackerDataImpl> StripID( DigitCollection ); 
    
    //Open status data
    LCCollectionVec * StatusCollection = dynamic_cast < LCCollectionVec * > (event->getCollection(_statusCollectionName)); 
    CellIDDecoder<TrackerRawDataImpl> StatusDecoder( StatusCollection ); 
    
    // We are assuming that the digit and status collections 
    // are aligned according to the sensorID. In other words, we are 
    // assuming the element i-th in the all the collections corresponds 
    // to the same detector. Let's test this ...

    streamlog_out( MESSAGE3 ) << "Found status collection " << endl;
    
    for ( size_t iDet = 0 ; iDet < DigitCollection->size(); ++iDet ) {

      streamlog_out( MESSAGE3 ) << "Found iDet " << iDet << endl;
       
      TrackerDataImpl * Digits = dynamic_cast<TrackerDataImpl* > ( DigitCollection->getElementAt(iDet) )  ;
      TrackerRawDataImpl * Status = dynamic_cast < TrackerRawDataImpl * >(StatusCollection->getElementAt(iDet));
      
      int dataID = static_cast<int> ( StripID(Digits)["sensorID"] );
      int statusID = static_cast<int> ( StatusDecoder(Status)["sensorID"] );
            
      if (dataID != statusID) {
        streamlog_out(ERROR3) << "Status collection not aligned to detector data!" << std::endl << std::endl;      
        //exit(-1); 
      }
      
    }
    
    _isStatusReady = true;
    
  } catch (  lcio::DataNotAvailableException ) {
    streamlog_out( MESSAGE2 ) << "Unable to to find status collection" << endl;
    _isStatusReady = false;
  }
   
}

} // Namespace



