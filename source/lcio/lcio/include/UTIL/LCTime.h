// -*- C++ -*-
#ifndef UTIL_LCTIME_H
#define UTIL_LCTIME_H 1

#include "LCIOTypes.h"

namespace UTIL {
  
  /** Helper class that allows to convert time stamps as defined in 
   *  LCEvent::getTimeStamp() ( ns since 1.1.1970 00:00:00 UTC) to
   *  ordinary date and time  and vice versa.<br>
   *  Dates and times are in UTC/GMT by default.<br>
   *  Negative time stamps are not supported and mapped to 0LL and 
   *  1.1.1970 00:00:00 UTC respectively.
   * 
   * @author gaede 
   * @version $Id: LCTime.h,v 1.4 2005/04/22 15:39:01 gaede Exp $
   */
  class LCTime {

  protected:
    /** Helper struct that holds the calendar time.
     */
    struct CalendarTime{
      short year  ;
      short  month ;
      short  day   ;
      short  hour ;
      short  min ;
      short  sec ;
      int ns ;
    } ;

  public: 
    
    /** Init with current time ( accuracy depends on OS/machine,etc.) in UTC/GMT.
     */
    LCTime() ;
    
    /** Init from 64 bit time stamp ( ns since 1.1.1970 00:00:00 ).
     */
    LCTime( EVENT::long64 time ) ;
    
    /** Init from date and time. 
     */
    LCTime(int year, int month, int day, int hour=0, int min=0, int s=0 ) ;
    
    /** Init from unix time ( s since 1.1.1970 00:00:00 ).
     */
    LCTime( int unixTime ) ;


    /** Unix time ( s since 1.1.1970 00:00:00 ).
     */
    int unixTime() { return ( _t / 1000000000LL ) ; } 

    /** 64bit time stamp as used in LCEvent::getTimestamp().
     */
    EVENT::long64 timeStamp() { return _t ; } 


    /** True if the year is a leap year in the Gregorian calendar:<br> 
     *  year is multiple of 4 and not multiple of 100<br>
     *  or year is multiple of 400.
     */
    bool isLeapYear( int year ) { return ( ( !( year % 4 ) &&  ( year%100 ) ) || !(year%400)   ) ; }

    /** The number if days in the given year in the Gregorian calendar. */
    int daysInYear( int year ) ;

    /** The number if days in the given month and year in the Gregorian calendar. */
    int daysInMonth( int month , int year ) ;

    /** Date in human readable format, e.g. 10.02.2005  10:54:29.123456789 :<br>
     *  dd.mm.yyyy  hh:mm:ss._ms_us_ns
     */
    std::string getDateString() ;

    /** Calendar time: year,month,day hour,min,sec and ns.
     */
    const CalendarTime& calendarTime() { return _d ; } 

    /** Calendar time:  year */
    int  year() { return  _d.year  ; }
    
    /** Calendar time:  month */
    int  month() { return _d.month   ; }
    
    /** Calendar time:  day */
    int  day() { return _d.day   ; }
    
    /** Calendar time:  hour */
    int  hour() { return _d.hour   ; }
    
    /** Calendar time:  min */
    int  min() { return _d.min   ; }
    
    /** Calendar time:  sec */
    int  sec() { return _d.sec  ; }
    
    /** Calendar time:  ns */
    int  ns() { return _d.ns   ; }
    
    /** Add the specified number of ns to the time.
     */
    void operator+=(EVENT::long64 t) ;


    virtual ~LCTime() { /*no_op*/; } 
    

    /** Tests the LCTime class with nDates random dates 
     *  Returns true if successful - throws Exception in case of error.
     */
    static bool test( int nDates ) ;
    
  protected: 
    
    EVENT::long64 _t ;  // time stamp in ns
    CalendarTime  _d ;  // calendar date and time
    
    
    // internal helper methods
    void convertToCalTime() ;

    void convertFromCalTime() ;

    static int dpm[13] ;  // days per month

  }; // class


} // namespace UTIL

#endif /* ifndef UTIL_LCTIME_H */
