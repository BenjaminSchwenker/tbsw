#ifndef UTIL_LCStdHepRdr_H
#define UTIL_LCStdHepRdr_H 1

#include "IMPL/LCCollectionVec.h"
#include "UTIL/lStdHep.hh"
#include "EVENT/LCIO.h"
#include "Exceptions.h"

namespace IMPL{

class LCEventImpl ;

}

namespace UTIL{
  
  /**Basic utility for reading a binary stdhep file and filling
   * a LCCollectionVec with MCParticles containing the stdhep
   * file information.
   * 
   * @author cassell
   * @version $Id: LCStdHepRdr.h,v 1.4 2007/11/12 16:39:04 gaede Exp $
   */
  class LCStdHepRdr{
    
  public:

	/** Open the stdhep input file in the constructer
	 */
    LCStdHepRdr(const char* evfile) ;

	/** noop
	 */
	~LCStdHepRdr() ;

    /** Read an event and return an LCCollectionVec of MCParticles.
     * @deprecated please use updateEvent()
     */
	IMPL::LCCollectionVec * readEvent() ;

    /** Reads the next stdhep event and adds a new MCParticle collection to the
     *  the event with default name 'MCParticle'
     * @throw IO::EndOfDataException if no event in stdhep file
     */
    void updateNextEvent( IMPL::LCEventImpl* evt , const char* colName=EVENT::LCIO::MCPARTICLE ) ;


    /** Print the file header to the given ostream.
     */
    void printHeader(std::ostream& os = std::cout ) ; 


    /** Return the charge of the particle times 3  - code copied from HepPDT package.
     */
    int threeCharge( int pdgID ) const ;

  private:
    
	lStdHep* _reader;
    

  }; // class

} // namespace UTIL

#endif /* ifndef UTIL_LCStdHepRdr_H */
