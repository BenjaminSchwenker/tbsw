#include "Time.hh"
#include "Exception.hh"
#include <ctime>
#include <iostream>


namespace eudaqinput {

  const std::string Time::DEFAULT_FORMAT = "%Y-%m-%d %H:%M:%S.%3";
  
  Time Time::Current() {
    timeval tv;
    if (gettimeofday(&tv, 0)) {
      EUDAQINPUT_THROWX(FileReadException, "Error getting current time" );
    }
    return tv;
  }



  Time::Time(int year, int month, int date, int hour, int minute, int sec, int usec) {
    struct tm time;
    //time.tm_gmtoff = 0;
    time.tm_isdst = -1;
    time.tm_year = year - 1900;
    time.tm_mon = month - 1;
    time.tm_mday = date;
    time.tm_hour = hour;
    time.tm_min = minute;
    time.tm_sec = sec + usec / 1000000;
    tv_sec = static_cast<long>(mktime(&time));
    tv_usec = usec % 1000000;
  }

  std::string Time::Formatted(const std::string & format) const {
    char buf[256];
    std::string fmt = format;
    size_t i, pos = 0;
    while ((i = fmt.find('%', pos)) != std::string::npos) {
      //std::cout << "i=" << i << std::endl;
      if (i >= fmt.length() - 1) break;
      pos = i+2;
      //std::cout << "pos=" << pos << std::endl;
      int c = fmt[i+1] - '0';
      if (c >= 1 && c <= 6) {
        //std::cout << "c=" << c << std::endl;
        std::string usecs = to_string(tv_usec, 6);
        //std::cout << "fmt=" << fmt << std::endl;
        fmt = std::string(fmt, 0, i) + std::string(usecs, 0, c) + std::string(fmt, i+2);
        //std::cout << "fmt=" << fmt << std::endl;
        pos += c - 2;
        //std::cout << "pos=" << pos << std::endl;
      }
    }
    time_t t = tv_sec;
    struct tm * tm = std::localtime(&t);
    std::strftime(buf, sizeof buf, fmt.c_str(), tm);
    return std::string(buf);
  }

}
