#ifndef Global_h
#define Global_h 1


namespace marlin{

  class StringParameters ;

  /** Simple global class for Marlin.
   *  Holds global parameters.
   *
   *  @author F. Gaede, DESY
   *  @version $Id: Global.h,v 1.4 2005-10-11 12:56:28 gaede Exp $ 
   * 
   *  B. Schwenker: removed gear dependency
   */
  class Global{
    
  public:
    
    static StringParameters* parameters ;
  };
  
  
} // end namespace marlin 
#endif
