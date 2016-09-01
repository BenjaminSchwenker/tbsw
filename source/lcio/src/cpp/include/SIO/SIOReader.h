#ifndef SIO_SIOREADER_H
#define SIO_SIOREADER_H 1

#include <string>
#include <set>
#include <map>
#include "IO/LCReader.h"
#include "IO/LCEventListener.h"
#include "IO/LCRunListener.h"

#include "IOIMPL/LCEventIOImpl.h"
#include "IOIMPL/LCRunHeaderIOImpl.h"
#include "LCIOTypes.h"

class SIO_record ;
class SIO_stream ;    


namespace SIO {

class SIORunHeaderHandler ;
class SIOEventHandler ;
  
/** Concrete implementation of LCWriter using SIO.
 * 
 * @author gaede
 * @version $Id: SIOReader.h,v 1.26 2008/12/10 08:10:59 gaede Exp $
 */
  class SIOReader : public IO::LCReader {
    
    typedef std::map< EVENT::long64 , EVENT::long64 > EventMap ;

  public:
    
    /** Default constructor.
     */
    SIOReader( int lcReaderFlag=0 ) ;
    
    // Destructor
    virtual ~SIOReader() ;
    

    /** Opens a list of files for reading (read-only). All subsequent
     * read operations will operate on the list, i.e. if an EOF is encountered
     * the next file in the list will be opened and read transparently to the
     * user.
     * @throws IOException
     */
    virtual void open(const std::vector<std::string>& filenames) 
      throw (IO::IOException, std::exception) ;


    /** Opens a file for reading (read-only).
     * @throws IOException
     */
    virtual void open(const std::string & filename) throw (IO::IOException, std::exception) ;
    
    /** Reads the next run header from the file. 
     *
     * @throws IOException
     */
    virtual EVENT::LCRunHeader * readNextRunHeader() throw (IO::IOException, std::exception) ;

    /** Same as readNextRunHeader() but allows to set the access mode 
     *  LCIO::READ_ONLY (default) or LCIO::Update. 
     *
     * @throws IOException
     */
    virtual EVENT::LCRunHeader * readNextRunHeader(int accessMode) throw (IO::IOException, std::exception) ;


    /** Reads the next event from the file. 
     *
     * @throws IOException
     */
    virtual EVENT::LCEvent* readNextEvent() throw (IO::IOException, std::exception) ;
    

    /** Same as readNextRunHeader() but allows to set the access mode 
     *  LCIO::READ_ONLY (default) or LCIO::Update
     *
     * @throws IOException
     */
    virtual EVENT::LCEvent* readNextEvent( int accessMode) throw (IO::IOException, std::exception) ;
    

    /** Skips the next n events from the current position. In fact simply reads the next n
     *  event headers so that the next event read is the (n+1)-th event.
     */
    virtual void skipNEvents(int n) ;


    /** Reads the specified event from file. 
     *  To be used with care: events have to be read in sequential 
     *  order (as LCIO has no direct access yet).
     *
     * @throws IOException
     */
    virtual EVENT::LCEvent * readEvent(int runNumber, int evtNumber) 
      throw (IO::IOException, std::exception/*, EVENT::NotAvailableException */) ;

    /** Closes the output file/stream etc.
     *
     * @throws IOException
     */
    virtual void close() throw (IO::IOException, std::exception) ;
    
    // interface for listeners
 
    /** Registers a listener for reading LCEvents from a stream.
     */ 
    virtual void registerLCEventListener(IO::LCEventListener * ls) ;

    /** Remove a listener for reading LCEvents from a stream.
     */ 
    virtual void removeLCEventListener(IO::LCEventListener * ls) ;

    /** Registers a listener for reading LCEventsLCRunHeaders from a stream.
     */ 
    virtual void registerLCRunListener(IO::LCRunListener * ls) ;

    /** Remove a listener for reading LCRunHeaders from a stream.
     */ 
    virtual void removeLCRunListener(IO::LCRunListener * ls) ;

    /** Reads the input stream and notifies registered 
     * listeners according to the object type 
     * found in the stream. 
     *
     * @throws IOException
     * @throws EndOfException
     */
    virtual void readStream() throw (IO::IOException, std::exception) ;

    /** Reads maxRecord from the input stream and notifies registered 
     * listeners according to the object type found in the stream. 
     * Throws EndOfException if less than maxRecord records are found in the stream. 
     *
     * @throws IOException
     * @throws EndOfException
     */
    virtual void readStream(int maxRecord) throw (IO::IOException, std::exception) ;




  protected:

    void setUpHandlers() ;
    void readRecord() throw (IO::IOException , IO::EndOfDataException , std::exception) ;


    void postProcessEvent() ;

    void getEventMap() ;

  protected:
    
    // we need an SIO record for every type
//     SIO_record *_evtRecord ;
//     SIO_record *_hdrRecord ;
//     SIO_record *_runRecord ;
    SIO_record *_dummyRecord ;  // used for reading arbitrary records
    SIO_stream *_stream ;

    SIORunHeaderHandler* _runHandler ;
    SIOEventHandler* _evtHandler ;

  private:
    
    IOIMPL::LCEventIOImpl *_defaultEvt ; // used to add collections when reading 
    IOIMPL::LCEventIOImpl **_evtP ;
    IOIMPL::LCRunHeaderIOImpl **_runP ;

    std::set<IO::LCRunListener*> _runListeners ;
    std::set<IO::LCEventListener*> _evtListeners ;
    
    const std::vector<std::string>* _myFilenames ;
    unsigned int _currentFileIndex ;

    EventMap _evtMap ;
    const bool _readEventMap ;

  }; // class
} // namespace

#endif /* ifndef SIO_SIOREADER_H */
