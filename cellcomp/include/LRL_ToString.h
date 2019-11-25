#ifndef LRL_TOSTRING_H
#define LRL_TOSTRING_H

#define DEBUG_REPORTS
#ifdef DEBUG_REPORTS
#define DEBUG_REPORT_STRING( s ) \
   std::cout << s << std::endl;
#else
#define DEBUG_REPORT_STRING( s )
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


//#pragma warning( disable : 4505) // unreferenced local function has been removed
/*
static std::string LRL_ToString( const char c ) {
   return std::string( 1, c );
}
*/

template<typename T>
std::string LRL_ToString( const T& t ) {
   std::ostringstream ostr;
   ostr << t;
   return ( ostr.str( ) );
}

template<typename T1, typename T2>
std::string LRL_ToString( const T1& t1, const T2& t2 ) {
   return LRL_ToString( t1 ) + " " + LRL_ToString( t2 );
}

template<typename T1, typename T2, typename T3>
std::string LRL_ToString( const T1& t1, const T2& t2, const T3& t3 ) {
   return LRL_ToString( t1 ) + " " + LRL_ToString( t2 ) + " " + LRL_ToString( t3 );
}

template<typename T1, typename T2, typename T3, typename T4>
std::string LRL_ToString( const T1& t1, const T2& t2, const T3& t3, const T4& t4 ) {
   return LRL_ToString( t1 ) + " " + LRL_ToString( t2 ) + " " + LRL_ToString( t3 ) + " "
      + LRL_ToString( t4 );
}

template<typename T1, typename T2, typename T3, typename T4, typename T5>
std::string LRL_ToString( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
   const T5& t5 ) {
   return LRL_ToString( t1 ) + " " + LRL_ToString( t2 ) + " " + LRL_ToString( t3 ) + " "
      + LRL_ToString( t4 ) + " " + LRL_ToString( t5 );
}

template<typename T1, typename T2, typename T3, typename T4, typename T5,
   typename T6>
   std::string LRL_ToString( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
   const T5& t5, const T6& t6 ) {
      return LRL_ToString( t1 ) + " " + LRL_ToString( t2 ) + " " + LRL_ToString( t3 ) + " "
         + LRL_ToString( t4 ) + " " + LRL_ToString( t5 ) + " " + LRL_ToString( t6 );
   }

template<typename T1, typename T2, typename T3, typename T4, typename T5,
   typename T6, typename T7>
   std::string LRL_ToString( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
   const T5& t5, const T6& t6, const T7& t7 ) {
      return LRL_ToString( t1 ) + " " + LRL_ToString( t2 ) + " " + LRL_ToString( t3 ) + " "
         + LRL_ToString( t4 ) + " " + LRL_ToString( t5 ) + " " + LRL_ToString( t6 ) + " "
         + LRL_ToString( t7 );
   }

template<typename T1, typename T2, typename T3, typename T4, typename T5,
   typename T6, typename T7, typename T8>
   std::string LRL_ToString( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
   const T5& t5, const T6& t6, const T7& t7, const T8& t8 ) {
      return LRL_ToString( t1 ) + " " + LRL_ToString( t2 ) + " " + LRL_ToString( t3 ) + " "
         + LRL_ToString( t4 ) + " " + LRL_ToString( t5 ) + " " + LRL_ToString( t6 ) + " "
         + LRL_ToString( t7 ) + " " + LRL_ToString( t8 );
   }

template<typename T1, typename T2, typename T3, typename T4, typename T5,
   typename T6, typename T7, typename T8, typename T9>
   std::string LRL_ToString(const T1& t1, const T2& t2, const T3& t3, const T4& t4,
      const T5& t5, const T6& t6, const T7& t7, const T8& t8, const T9& t9) {
   return LRL_ToString(t1) + " " + LRL_ToString(t2) + " " + LRL_ToString(t3) + " "
      + LRL_ToString(t4) + " " + LRL_ToString(t5) + " " + LRL_ToString(t6) + " "
      + LRL_ToString(t7) + " " + LRL_ToString(t8) + " " + LRL_ToString(t9);
}

template<typename T1, typename T2, typename T3, typename T4, typename T5,
   typename T6, typename T7, typename T8, typename T9, typename T10>
   std::string LRL_ToString(const T1& t1, const T2& t2, const T3& t3, const T4& t4,
      const T5& t5, const T6& t6, const T7& t7, const T8& t8, const T9& t9, const T10& t10) {
   return LRL_ToString(t1) + " " + LRL_ToString(t2) + " " + LRL_ToString(t3) + " "
      + LRL_ToString(t4) + " " + LRL_ToString(t5) + " " + LRL_ToString(t6) + " "
      + LRL_ToString(t7) + " " + LRL_ToString(t8) + " " + LRL_ToString(t9)
      + " " + LRL_ToString(t10);
}

template<typename T>
static std::string LRL_ToString( const std::vector<T>& t ) {
   std::string s;
   const size_t n = (size_t)t.size( );
   for( size_t i = 0; i < n; ++i ) {
      s += LRL_ToString(t[i]) + "\t";
   }
   return s;
}

template<typename T>
std::string LRL_ToString( const std::string& s, const T& t ) {
   return ( s + LRL_ToString( t ) );
}

template<typename T>
std::string LRL_ToString( const T& t, const std::string& s ) {
   return ( LRL_ToString( t ) + s );
}

template<typename T>
std::string ToQuotedString( const T& t ) {
   return ("\"" + LRL_ToString(t) +"\"" );
}

#endif  // LRL_TOSTRING_H

#undef MAKE_TOSTRING

