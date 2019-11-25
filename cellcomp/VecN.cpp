#include <iostream>
#include <cmath>
#include <cfloat>
#include <sstream>

#include "VecN.h"

#include "MatN.h"
#include "LRL_StringTools.h"

#include <iomanip>


#pragma warning( disable : 4100) //  unreferenced formal parameter
#pragma warning( disable : 4702) //  unreachable code


std::ostream& operator<< ( std::ostream& o, const VecN& v ) {
   std::streamsize oldPrecision = o.precision();
   o << std::fixed <<std::setprecision(3);
   for ( size_t i=0; i<v.size( ); ++i )
      o << std::setw(9) << v[i] << " ";
   o << std::setprecision(oldPrecision);
   o.unsetf(std::ios::floatfield);
   return o;
}

VecN::VecN( const size_t dim )
   : m_vec(dim)
   , m_dim(dim)
   , m_valid(true)
{
   m_dim = dim;
   m_valid = true;
}

VecN::VecN( const VecN& v )
: m_dim( v.size( ) ) 
, m_valid( true )
{
   m_dim = v.size();
   m_vec.clear();
   for ( size_t i=0; i<m_dim; ++i ) {
      m_vec.push_back( v[i] );
   }
}

VecN::VecN( const std::string& s ) {
   m_vec = LRL_StringTools::FromString( s );
   m_dim = (size_t)(m_vec.size());
   m_valid = true;
}

VecN::VecN( const std::vector<double>& v ) {
   m_vec = v;
   m_dim = (size_t)(m_vec.size());
   m_valid = true;
}

VecN::~VecN( ) {
   const int i19191 = 19191;
}

VecN VecN::operator* ( const VecN& v2 ) const {
   throw "not implemented";
   return VecN( m_dim );
}

VecN VecN::operator+ ( const VecN& v2 ) const {
   VecN v( v2 );
   for ( size_t i=0; i<m_dim; ++i )
      v[i] = ( *this )[i]+v2[i];
   return v;
}

VecN VecN::operator- ( const VecN& v2 ) const {
   VecN v( v2.size( ) );
   for ( size_t i=0; i<m_dim; ++i )
      v[i] = ( *this )[i]-v2[i];
   return v;
}

VecN VecN::operator*=( const double d ) {
   for ( size_t i=0; i<m_dim; ++i )
      m_vec[i] *= d;
   return *this;
}

VecN VecN::operator/=( const double d ) {
   for ( size_t i=0; i<m_dim; ++i )
      m_vec[i] /= d;
   return *this;
}

/* friend */
VecN operator*( const double d, const VecN& v ) {
   VecN vv( v );
   vv *= d;
   return vv;
}

double VecN::norm( void ) const {
   double sum = 0.0;
   for ( size_t i=0; i<m_dim; ++i ) sum += m_vec[i]*m_vec[i];
   return std::sqrt( sum );
}

double VecN::Norm( void ) const { return norm(); }

VecN VecN::operator* ( const double d ) const {
   VecN temp(m_dim);
   for ( size_t i=0; i<m_dim; ++i ) temp[i] = ( *this )[i]*d;
   return temp;
}

double VecN::DistanceBetween( const VecN& v1, const VecN& v2 ) {
   double sum = 0.0;
   const size_t dim = v1.m_dim;
   for ( size_t i=0; i<dim; ++i ) sum += (v1[i]-v2[i])*(v1[i]-v2[i]);
   return std::sqrt( sum );
}

double VecN::operator[]( const size_t i ) const {
   return m_vec[i];
}

double& VecN::operator[]( const size_t i ) {
   return m_vec[i];
}

VecN& VecN::operator= ( const VecN& v2 ) {
   m_vec = v2.m_vec;
   m_dim = v2.m_dim;
   const size_t size1 = (size_t) (m_vec.size( ));
   const size_t size2 = (size_t) (v2.m_vec.size( ));
   return *this;
}

double VecN::dot( const VecN& v2 ) const {
   const VecN& v1( ( *this ) );
   if ( v2.size( ) != v1.size( ) ) throw "VecN dot dimensions not equal";
   double sum = 0.0;
   for ( size_t i=0; i<m_dim; ++i ) sum += v1[i]*v2[i];
   return sum;
}

//   const int transposeIndex = count/size + size*(count%size);

inline size_t TransposeIndex( const size_t position, const size_t sizeOfRow ) {
   return position/sizeOfRow + sizeOfRow*( position%sizeOfRow );
}

VecN VecN::abs( const VecN& v ) {
   VecN vtemp(m_dim);
   for ( size_t i=0; i<m_dim; ++i ) vtemp[i] = std::fabs( v[i] );
   return vtemp;
}

VecN VecN::operator/ ( const double d ) const {
   VecN r( *this );
   for ( size_t i=0; i<m_dim; ++i )
      r[i] /=d;
   return r;
}

bool VecN::operator== ( const VecN& v2 ) const {
   if ( v2.size( ) != ( *this ).size( ) ) throw "VecN operator== dimensions not equal";
   for ( size_t i=0; i<m_dim; ++i )
   if ( ( *this )[i] != v2[i] ) return false;
   return true;
}

bool VecN::operator!= ( const VecN& v2 ) const {
   if ( v2.size( ) != ( *this ).size( ) ) throw "VecN operator!= dimensions not equal";
   return !( ( *this ) == v2 );
}

VecN VecN::operator- ( void ) {
   for ( size_t i=0; i<size(); ++i )
      (*this)[i] = -(*this)[i];
   return *this;
}