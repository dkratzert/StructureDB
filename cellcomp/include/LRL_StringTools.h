#ifndef LRL_StringTools_H
#define LRL_StringTools_H
#include <string>
#include <vector>

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <iterator>
#include <sstream>

class LRL_StringTools
{
public:

   static std::string strToupper(const std::string& s);
   static std::vector<double> FromString(const std::string& s);


};



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
std::string LRL_StringTools::strToupper(const std::string& s) {
   std::string ss(s);
   for (size_t i = 0; i < ss.length(); ++i)
	   ss[i] = static_cast<char>(toupper(ss[i]));
   return(ss);
}

std::vector<double> LRL_StringTools::FromString(const std::string& s) {
   std::istringstream istr(s);
   std::vector<double> t;
   double d;
   int i = 0;
   while (istr && !istr.eof()) {
      istr >> d;
      if (istr.eof() && d == DBL_MAX) break;
      t.push_back(d);
      ++i;
      d = DBL_MAX;
   }

   return t;
}



#endif  // LRL_StringTools_H
