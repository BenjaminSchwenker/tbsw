#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>

#include "marlin/DataSourceProcessor.h"
#include "marlin/EventModifier.h"
#include "streamlog/streamlog.h"
#include "streamlog/logbuffer.h"

#include <time.h>

namespace marlin{

ProcessorMgr* ProcessorMgr::_me = 0 ;

static clock_t start_t , end_t ;
typedef std::map< Processor* , std::pair< double  , int > > TimeMap ;
static TimeMap tMap ;



// helper for sorting procs wrt to processing time
struct Cmp{
    bool operator()(const TimeMap::value_type& v1, const TimeMap::value_type& v2 ) {
        // inverse sort:
        return v1.second.first > v2.second.first ;
    }
} ;

struct ProcMgrStopProcessing : public StopProcessingException {
    ProcMgrStopProcessing(const std::string m){
        StopProcessingException::message = m  ;
    }
};



// create a dummy streamlog stream for std::cout
streamlog::logstream my_cout ;


ProcessorMgr* ProcessorMgr::instance() {

    if( _me == 0 )
        _me = new ProcessorMgr ;

    return _me ;
}

ProcessorMgr::ProcessorMgr(){
    _time_outsideProcessors=0;
    _time_startup=0;
    _time_end=0;
    _isFirstEvent=true;
    end_t=clock();
}



void ProcessorMgr::registerProcessor( Processor* processor ){

    const std::string& name = processor->type()  ;

    if( _map.find( name ) != _map.end() ){

        //     std::cerr << " ProcessorMgr::registerProcessor: processor " <<  name
        // 	      << " already registered ! "
        // 	      << std::endl   ;

        return ;
    }
    else

        _map[ name ] = processor ;

}

void ProcessorMgr::readDataSource( int numEvents ) {

    for(  ProcessorList::iterator it = _list.begin() ;
          it != _list.end() ; it++ ){

        DataSourceProcessor* dSP = dynamic_cast<DataSourceProcessor*>( *it ) ;

        if( dSP != 0 )
            dSP->readDataSource( numEvents ) ;

    }
}


void  ProcessorMgr::dumpRegisteredProcessors() {

    typedef ProcessorMap::iterator MI ;
    
    std::cout  << "  ##########################################" << std::endl
               << "  #                                        #" << std::endl
               << "  #     Example steering file for marlin   #" << std::endl
               << "  #                                        #" << std::endl
               << "  ##########################################" << std::endl
               <<  std::endl ;
    
    std::cout  << ".begin Global  ---------------------------------------" << std::endl
               << "   LCIOInputFiles simjob.slcio " << std::endl
               << std::endl
               << "  # the active processors that are called in the given order" << std::endl
               << "   ActiveProcessors MyAIDAProcessor" << std::endl
               << "   ActiveProcessors MyTestProcessor" << std::endl
               << "   ActiveProcessors MyLCIOOutputProcessor" << std::endl
               << std::endl
               << "  # limit the number of processed records (run+evt):" << std::endl
               << "   MaxRecordNumber 5001" << std::endl
               << std::endl
               << "  # skip the first  n events  " << std::endl
               << "  SkipNEvents  0 " << std::endl
               << "  # don't call the check method of the processors if \"true\"" << std::endl
               << "   SupressCheck false" << std::endl
               << ".end   -----------------------------------------------" << std::endl
               <<  std::endl
                <<  std::endl ;
    
    
    for(MI i=_map.begin() ; i!= _map.end() ; i++) {
        i->second->printDescription() ;
    }
}
void  ProcessorMgr::dumpRegisteredProcessorsXML() {

    typedef ProcessorMap::iterator MI ;
    
    std::cout  <<  std::endl
                << "<marlin xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                << "xsi:noNamespaceSchemaLocation=\"http://ilcsoft.desy.de/marlin/marlin.xsd\">"
                <<  std::endl ;
    
    std::cout  <<  " <global>" << std::endl
               <<  "  <parameter name=\"LCIOInputFiles\" value=\"\" />" << std::endl
               <<  "  <parameter name=\"MaxRecordNumber\" value=\"-1\" />  " << std::endl
               <<  "  <parameter name=\"SkipNEvents\" value=\"0\" />  " << std::endl
               <<  "  <parameter name=\"SupressCheck\" value=\"false\" />  " << std::endl
               <<  "  <parameter name=\"GearXMLFile\" value=\"\" />  " << std::endl
               <<  "  <parameter name=\"Verbosity\" options=\"DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT\" value=\"MESSAGE3\" /> " << std::endl
               <<  " </global>" << std::endl
               << std::endl ;

    for(MI i=_map.begin() ; i!= _map.end() ; i++) {
        i->second->printDescriptionXML() ;
    }

    std::cout  <<  std::endl
                << "</marlin>"
                <<  std::endl ;

}

std::set< std::string > ProcessorMgr::getAvailableProcessorTypes(){

    std::set< std::string > ptypes;

    for(ProcessorMap::iterator i=_map.begin() ; i!= _map.end() ; i++) {
        ptypes.insert(i->first);
    }
    return ptypes;
}

Processor* ProcessorMgr::getProcessor( const std::string& type ){
    return _map[ type ] ;
}

Processor* ProcessorMgr::getActiveProcessor( const std::string& name ) {
    return _activeMap[ name ] ;
}

void ProcessorMgr::removeActiveProcessor(  const std::string& name ) {


    _list.remove( _activeMap[name] ) ;
    _activeMap.erase( name ) ;

}


bool ProcessorMgr::addActiveProcessor( const std::string& processorType ,
                                       const std::string& processorName ,
                                       StringParameters* parameters ,
                                       const std::string condition) {

    Processor* processor = getProcessor( processorType ) ;



    if( processor == 0 ) {
        std::stringstream sstr ;
        sstr << " ProcessorMgr::registerProcessor: unknown processor with type " <<  processorType  << " ! " << std::endl   ;
        throw Exception( sstr.str() );
    }

    
    if( _activeMap.find( processorName ) != _activeMap.end() ){

        std::cerr << " ProcessorMgr::addActiveProcessor: processor " <<  processorName
                  << " already registered ! "
                  << std::endl ;
        return false ;

    } else {

        Processor* newProcessor = processor->newProcessor() ;
        newProcessor->setName( processorName ) ;
        _activeMap[ processorName ] = newProcessor ;
        _list.push_back( newProcessor ) ;
        _conditions.addCondition( processorName, condition ) ;

        if( parameters != 0 ){
            newProcessor->setParameters( parameters  ) ;
        }
        //       // keep a copy of the output processor
        //       if( processorType == "LCIOOutputProcessor" ){
        // 	_outputProcessor = dynamic_cast<LCIOOutputProcessor*>( newProcessor ) ;
        //       }
    }

    return true ;
}


void ProcessorMgr::init(){

    streamlog::logbuffer* lb = new streamlog::logbuffer( std::cout.rdbuf() ,  &my_cout ) ;
    std::cout.rdbuf(  lb ) ;

    //     for_each( _list.begin() , _list.end() , std::mem_fun( &Processor::baseInit ) ) ;
    
    for( ProcessorList::iterator it = _list.begin() ; it != _list.end() ; ++it ) {

        streamlog::logscope scope( streamlog::out ) ; scope.setName(  (*it)->name()  ) ;
        streamlog::logscope scope1(  my_cout ) ; scope1.setName(  (*it)->name()  ) ;

        (*it)->baseInit() ;

        tMap[ *it ] = std::make_pair( 0 , 0 )  ;
    }

}

void ProcessorMgr::processRunHeader( LCRunHeader* run){
    
    //     for_each( _list.begin() , _list.end() ,  std::bind2nd(  std::mem_fun( &Processor::processRunHeader ) , run ) ) ;
    for( ProcessorList::iterator it = _list.begin() ; it != _list.end() ; ++it ) {

        streamlog::logscope scope( streamlog::out ) ; scope.setName(  (*it)->name()  ) ;
        streamlog::logscope scope1(  my_cout ) ; scope1.setName(  (*it)->name()  ) ;

        (*it)->processRunHeader( run ) ;
    }
}


//   void ProcessorMgr::modifyEvent( LCEvent * evt ) { 
//     if( _outputProcessor != 0 )
//       _outputProcessor->dropCollections( evt ) ;
//   }



void ProcessorMgr::modifyEvent( LCEvent* evt ){

    static bool first = true  ;
    typedef std::vector<EventModifier*> EMVec ;
    static  EMVec  emv ;

    if( first ) {

        for(  ProcessorList::iterator it = _list.begin() ;
              it != _list.end() ; it++ ){

            EventModifier* em = dynamic_cast<EventModifier*>( *it ) ;

            if( em != 0 ) {
                emv.push_back( em ) ;
                streamlog_out( WARNING4 ) << " -----------   " << std::endl
                                          << " the following processor will modify the LCIO event :  "
                                          << (*it)->name()  << " !! " <<  std::endl
                                          << " ------------  "   << std::endl ;
            }
        }

        first = false ;
    }

    for( EMVec::iterator it = emv.begin();  it !=  emv.end()  ; ++ it) {

        streamlog::logscope scope( streamlog::out ) ; scope.setName(  (*it)->name()  ) ;

        streamlog::logscope scope1(  my_cout ) ; scope1.setName(  (*it)->name()  ) ;

        (*it)->modifyEvent( evt ) ;
    }
    //    for_each( emv.begin() , emv.end() ,   std::bind2nd(  std::mem_fun( &EventModifier::modifyEvent ) , evt ) ) ;

}


void ProcessorMgr::processEvent( LCEvent* evt ){

    //     static bool isFirstEvent = true ;

    //     for_each( _list.begin() , _list.end() ,   std::bind2nd(  std::mem_fun( &Processor::processEvent ) , evt ) ) ;

    //     if( Global::parameters->getStringVal("SupressCheck") != "true" ) {

    //       for_each( _list.begin() , _list.end(),
    // 		std::bind2nd( std::mem_fun( &Processor::check ) , evt ) ) ;
    //     }
    
    //     if ( isFirstEvent ) {
    //       isFirstEvent = false;
    //       for_each( _list.begin(), _list.end() ,
    // 		std::bind2nd( std::mem_fun( &Processor::setFirstEvent ),isFirstEvent )) ;
    //     }

    _conditions.clear() ;

    bool check = ( Global::parameters->getStringVal("SupressCheck") != "true" ) ;

    if(_isFirstEvent) _time_startup=double(clock()-end_t);
    try{

        for( ProcessorList::iterator it = _list.begin() ; it != _list.end() ; ++it ) {

            if( _conditions.conditionIsTrue( (*it)->name() ) ) {

                streamlog::logscope scope( streamlog::out ) ; scope.setName(  (*it)->name()  ) ;
                streamlog::logscope scope1(  my_cout ) ; scope1.setName(  (*it)->name()  ) ;

                start_t =  clock () ;  // start timer
                _time_outsideProcessors+=double(start_t-end_t);

                (*it)->processEvent( evt ) ;
                if( check )  (*it)->check( evt ) ;
                end_t =  clock () ;  // stop timer


                TimeMap::iterator itT = tMap.find( *it ) ;

                itT->second.first += double( end_t - start_t  ) ;
                itT->second.second ++ ;


                (*it)->setFirstEvent( false ) ;
            }
        }
    } catch( SkipEventException& e){

        ++ _skipMap[ e.what() ] ;
    }
}


void ProcessorMgr::setProcessorReturnValue( Processor* proc, bool val ) {

    _conditions.setValue( proc->name() , val ) ;

}
void ProcessorMgr::setProcessorReturnValue( Processor* proc, bool val,
                                            const std::string& name){
    
    std::string valName = proc->name() + "." + name ;
    _conditions.setValue( valName , val ) ;
}

void ProcessorMgr::end(){

    //     for_each( _list.begin() , _list.end() ,  std::mem_fun( &Processor::end ) ) ;

    //    for_each( _list.rbegin() , _list.rend() ,  std::mem_fun( &Processor::end ) ) ;

    for( ProcessorList::reverse_iterator it = _list.rbegin() ; it != _list.rend() ; ++it ) {

        streamlog::logscope scope( streamlog::out ) ; scope.setName(  (*it)->name()  ) ;
        streamlog::logscope scope1(  my_cout ) ; scope1.setName(  (*it)->name()  ) ;

        (*it)->end() ;
    }
    //     if( _skipMap.size() > 0 ) {
    streamlog_out(MESSAGE3)  << " --------------------------------------------------------- " << std::endl
                             << "  Events skipped by processors : " << std::endl ;
    
    unsigned nSkipped = 0 ;
    for( SkippedEventMap::iterator it = _skipMap.begin() ; it != _skipMap.end() ; it++) {

        streamlog_out(MESSAGE3) << "       " << it->first << ": \t" <<  it->second << std::endl ;

        nSkipped += it->second ;
    }
    streamlog_out(MESSAGE3)  << "  Total: " << nSkipped  << std::endl ;
    streamlog_out(MESSAGE3)  << " --------------------------------------------------------- "
                             << std::endl
                             << std::endl ;
    //     }
    
    // ----- print timing information ----------
    
    streamlog_out(MESSAGE3)  << " --------------------------------------------------------- " << std::endl
                             << "      Time used by processors ( in processEvent() ) :      " << std::endl
                             << std::endl ;
    

    
    

    //    for( ProcessorList::iterator it = _list.begin() ; it != _list.end() ; ++it ) {
    //      TimeMap::iterator itT = tMap.find( *it ) ;

    // do not sort procs wrt processing time :
    typedef std::list< TimeMap::value_type > TMList  ;
    TMList l ;
    std::copy(  tMap.begin() , tMap.end() , std::back_inserter( l ) )  ;
    //l.sort( Cmp() ) ;
    
    double tTotal = 0.0 ;
    int evtTotal = 0 ;

    for( TMList::iterator itT = l.begin() ; itT != l.end() ; ++ itT ) {

        // fg: this does not work  w/ streamlog !?
        //      streamlog_out(MESSAGE) << std::ios_base::left << std::setw(30)  <<  (*it)->name() ;
        // copy string to fixed size char* ->

        char cName[40] = "                           ";
        const std::string& sName = itT->first->name()  ;
        unsigned nChar = ( sName.size() > 25 ?  25 : sName.size() )  ;
        for(unsigned  i=0 ; i< nChar ; i++ ) {
            cName[i] = sName[i] ;
        }


        double tProc = itT->second.first  / double(CLOCKS_PER_SEC) ;

        tTotal += tProc ;

        int evtProc = itT->second.second ;

        if( evtProc > evtTotal )
            evtTotal = evtProc ;

        streamlog_out(MESSAGE3)  <<  cName
                                  <<  std::setw(10)   <<std::setprecision (3)  <<tProc  << " s in "
                                   <<  std::setw(10) << evtProc << " events  ==> " ;

        if( evtProc > 0 ){
            streamlog_out(MESSAGE3)  <<  std::setw(10) <<std::setprecision (6)   << 1000.*tProc / evtProc << " [ ms/evt.] "  ;
        }else{
            streamlog_out(MESSAGE3)  <<  std::setw(10)   << "NaN"  << " [ s/evt.] "  ;
        }
        streamlog_out(MESSAGE3)  <<  std::endl ;

    }
    _time_outsideProcessors=_time_outsideProcessors/ double(CLOCKS_PER_SEC);
    _time_startup=_time_startup/ double(CLOCKS_PER_SEC);
    tTotal+=_time_outsideProcessors;
    streamlog_out(MESSAGE3)  <<  "Time spent outside processor         "
                              <<  std::setw(10)    <<std::setprecision (3) << _time_outsideProcessors << " s in "
                               <<  std::setw(10) << evtTotal << " events  ==> " ;
    if( evtTotal > 0 ){
        streamlog_out(MESSAGE3)  <<  std::setw(10) <<std::setprecision (6) << 1000.*_time_outsideProcessors / evtTotal << " [ ms/evt.] "  ;
    }else{
        streamlog_out(MESSAGE3)  <<  std::setw(10)   << "NaN"  << " [ s/evt.] "  ;
    }
    streamlog_out(MESSAGE3)  <<  std::endl ;
    streamlog_out(MESSAGE3)  <<  "Time for startup                     "
                              <<  std::setw(10)    <<std::setprecision (3) << _time_startup <<  std::endl ;

    _time_end=double(clock()-end_t)/ double(CLOCKS_PER_SEC);
    streamlog_out(MESSAGE3)  <<  "Time spent in end()                  "
                              <<  std::setw(10)   <<std::setprecision (3)  << _time_end << std::endl;


    streamlog_out(MESSAGE3)  <<  "            Total:                   "
                              <<  std::setw(10)   <<std::setprecision (3)  << tTotal << " s in "
                               <<  std::setw(10) << evtTotal << " events  ==> " ;
    if( evtTotal > 0 ){
        streamlog_out(MESSAGE3)  <<  std::setw(10) <<std::setprecision (6) << 1000.*tTotal / evtTotal << " [ ms/evt.] "  ;
    }else{
        streamlog_out(MESSAGE3)  <<  std::setw(10)   << "NaN"  << " [ s/evt.] "  ;
    }
    streamlog_out(MESSAGE3)  <<  std::endl ;
    
    
    streamlog_out(MESSAGE3) << " --------------------------------------------------------- "  << std::endl ;
    
}


} // namespace marlin
