


#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <cmath>
#include <sstream>

#include "MatN.h"

MatN::MatN( const size_t dim/*=Default36*/ )
: m_dim( dim )
, m_rowdim( int( sqrt( double( dim ) ) ) ) {
   for ( size_t i=0; i<dim; ++i )
      m_mat.push_back( 19191.0 );
}

MatN::MatN( const MatN& m )
: m_dim( m.m_dim )
, m_mat( m.m_dim )
, m_rowdim( int( sqrt( double( m_dim ) ) ) ) {
   for ( size_t i=0; i<m_dim; ++i ) m_mat[i] = m.m_mat[i];
}

MatN::MatN( const std::string& s ) {
   std::istringstream istr( s );
   double d;
   std::vector<double> v;
   while ( istr && !istr.eof( ) ) {
      istr >> d;
      v.push_back( d );
   }

   m_mat = v;
   m_dim = (size_t)(m_mat.size( ));
   m_rowdim = int( sqrt( double( m_dim ) ) );

   if ( m_dim != m_rowdim*m_rowdim ) throw "bad dimension in constructor";
   const size_t rowsize = (size_t)std::sqrt( double(m_dim) );
   if ( m_mat.size( ) != rowsize*rowsize ) throw "string contains different number of element than dimension";
}

MatN MatN::operator* ( const MatN& m2 ) const {
   MatN m( ( *this ).GetDim( ) );
   if ( m.size( ) != m2.size( ) ) throw "bad dimensions in operator*";
   m.zeros( );
   const size_t size = (size_t)std::sqrt( double(m_dim) );
   const MatN& m1( *this );
   for ( size_t count=0; count<m_dim; ++count ) {
      for ( size_t i=0; i<size; ++i ) {
         m[count] += m1[size*( count/size )+i] *m2[( count%size )+size*i];
      }
   }
   return m;
}

MatN MatN::operator*= ( const MatN& m2 ) {
   if ( size( ) != m2.size( ) ) throw "bad dimensions in operator*=";
   const size_t size = (size_t)std::sqrt( double(m_dim) );
   const MatN m1( *this );
   (*this).zeros( );
   for ( size_t count=0; count<m_dim; ++count ) {
      for ( size_t i=0; i<size; ++i ) {
         (*this)[count] += m1[size*( count/size )+i] *m2[( count%size )+size*i];
      }
   }
   return *this;
}

MatN MatN::operator+ ( const MatN& m2 ) const {
   MatN m( ( *this ).size( ) );
   if ( m.size( ) != m2.size( ) ) throw "bad dimensions in operator+";
   for ( size_t i=0; i<m_dim; ++i )
      m[i] = ( *this )[i]+m2[i];
   return m;
}

MatN MatN::operator- (const MatN& m2) const {
   MatN m(*this);
   if (m.size() != m2.size()) throw "bad dimensions in operator-";
   for (size_t i = 0; i<m_dim; ++i)
      m[i] -= m2[i];
   return m;
}

MatN MatN::operator- (void) const { // unary
   MatN m(*this);
   for (size_t i = 0; i<m_dim; ++i)
      m[i] = -m[i];
   return m;
}

MatN MatN::transpose( void ) const {
   //  transpose a symmetrical matrix
   MatN m( m_dim );
   const size_t rowsize = ( *this ).m_rowdim;
   for ( size_t count=0; count<m_dim; ++count ) {
      const size_t transposeIndex = count/rowsize + rowsize*( count%rowsize );
      if ( count >= transposeIndex ) {
         m[transposeIndex] = ( *this )[count];
         m[count] = ( *this )[transposeIndex];
      }
   }
   return m;
}

VecN MatN::operator* ( const VecN& v ) const {
   VecN vout( v.size( ) );
   if ( ( *this ).GetRowDim( ) != v.GetDim( ) ) throw "bad dimensions in operator+";
   const size_t size = (size_t)std::sqrt( double(m_dim) );
   for ( size_t i=0; i<size; ++i ) vout[i] =0.0;

   int count = 0;
   for ( size_t im=0; im<m_dim; im+=size ) {
      double d = 0.0;
      for ( size_t j=0; j<size; ++j ) {
         d += ( *this )[im+j]* v[j];
      }
      vout[count] = d ;
      ++count;
   }

   return vout;
}

double MatN::norm ( void ) const {
   double sum = 0.0;
   for ( size_t i=0; i<m_dim; ++i )
      sum += m_mat[i]*m_mat[i];
   return std::sqrt( sum );
}

MatN MatN::operator* ( const double d ) const {
   MatN m2( ( *this ).size( ) );
   for ( size_t i=0; i<m_dim; ++i )
      m2.m_mat[i] = d*m_mat[i];
   return m2;
}

MatN MatN::operator/ ( const double d ) const {
   MatN m2( ( *this ).size( ) );
   for ( size_t i=0; i<m_dim; ++i )
      m2.m_mat[i] = m_mat[i]/d;
   return m2;
}

MatN operator* ( const double d, const MatN& m ) {
   return m*d;
}

MatN operator/ ( const double d, const MatN& m ) {
   return m/d;
}

MatN MatN::Eye( const MatN& m_in ) const {
   MatN m(m_in);
   const size_t rowsize = (size_t)( std::sqrt( double(m.size( ) ) ) );
   for ( size_t i=0; i<m_dim; ++i )
      m.m_mat[i] = 0.0;
   for ( size_t i=0; i<m_dim; i+=+rowsize+1 )
      m.m_mat[i] = 1.0;
   return m;
}

MatN MatN::Eye( const size_t n ) {
   MatN m(n);
   return m.Eye(m);
}

void MatN::Eye(void) {
   (*this) = Eye(*this);
}


bool MatN::IsUnit( void ) const {
   const MatN& m(*this);
   for ( size_t i=0; i<m_dim; i+=m_rowdim+1) 
      if ( m[i] !=1.0 ) return false;
   double sum = 0.0;
   for ( size_t i=0; i<m_dim; ++i ) sum += std::fabs(m[i]);
   return sum == m_rowdim;
}

double MatN::operator[]( const size_t i ) const {
   return m_mat[i];
}

double& MatN::operator[]( const size_t i ) {
   return m_mat[i];
}

MatN& MatN::operator= ( const MatN& m2 ) {
   for ( size_t i=0; i<m_dim; ++i ) m_mat[i] = m2.m_mat[i];
   return *this;
}

void MatN::zeros( ) {
   for ( size_t i=0; i<m_dim; ++i ) m_mat[i] = 0.0;
}

bool MatN::operator== ( const MatN& m2 ) const {
   bool b = true;
   for ( size_t i=0; i<m_dim; ++i ) b = b && m2[i]==( *this )[i];
   return b;
}

MatN MatN::FromString( const std::string& s ) {
   std::istringstream istr( s );
   MatN t( 1000 );
   double d;
   int i = 0;
   while ( istr && ! istr.eof() ) {
      istr >> d;
      t[i] = d;
      ++i;
   }
   (*this).SetDim(i);
   SetRowDim((size_t)(sqrt(i)));
   t.m_dim = i;
   t.SetRowDim((*this).m_rowdim);
   t.m_mat.resize(m_dim);
   return t;
}

inline size_t LinearIndex( const size_t row, const size_t col, const size_t sizeOfRow ) {
   return row*sizeOfRow + col;
}

std::ostream& operator<< ( std::ostream& o, const MatN& m ) {
   const size_t dim = m.GetDim( );
   for ( size_t i=0; i<dim; ++i ) {
      if ( i > 0 && i%m.GetRowDim( )==0 ) o << std::endl;
      o << ( std::fabs( m[i] ) < 1.0E-20? 0.0: m[i] ) << " ";
   }
   return o;
}

MatN Eye( const MatN& m ) { // this is a free, local function only
   throw("this no longer works and needs to be fleshed out explicitly");
   return Eye( m );
}

double MatN::operator() (const size_t i) {
   return m_mat[i];
}

double MatN::operator() (const size_t row, const size_t col) {
   return m_mat[LinearIndex(row, col, m_rowdim)];
}
