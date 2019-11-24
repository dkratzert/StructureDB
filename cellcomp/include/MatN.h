#ifndef MAT_N_H
#define MAT_N_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "VecN.h"

class VecN;
class MatN;

MatN operator* ( const double d, const MatN& m );
MatN operator/ ( const double d, const MatN& m );
MatN Eye( const MatN& m );

class MatN {
public:
   friend std::ostream& operator<< ( std::ostream&, const MatN& );

   MatN( const size_t dim );
   MatN( const std::string& s );
   MatN( const MatN& m );
   virtual ~MatN( void ) {};

   virtual MatN operator* ( const MatN& m2 ) const;
   virtual MatN operator*= ( const MatN& m2 );
   virtual MatN operator+ ( const MatN& m2 ) const;
   virtual MatN operator- (const MatN& m2) const;
   virtual MatN operator- (void) const; // unary
   virtual MatN transpose( void ) const;

   virtual size_t size( void ) const { return m_dim; }

   virtual VecN operator* ( const VecN& v ) const;
   virtual double norm( void ) const;

   virtual MatN operator* ( const double d ) const;
   virtual MatN operator/ ( const double d ) const;
   virtual double operator[]( const size_t i ) const;
   virtual double& operator[]( const size_t i );
   virtual MatN& operator= ( const MatN& m2 );
   virtual bool operator== ( const MatN& m2 ) const;
   virtual void zeros( void );
   virtual MatN FromString( const std::string& s );
   virtual size_t GetDim(void) const { return m_dim; }
   virtual size_t GetRowDim(void) const { return m_rowdim; }
   virtual size_t GetColDim(void) const { return m_rowdim; }
   virtual void SetDim(const size_t n) { m_dim = n; }
   virtual void SetRowDim(const size_t n) { m_rowdim = n; }
   virtual void SetColDim(const size_t n) { m_rowdim = n; }

   virtual MatN Eye( const MatN& m_in ) const;
   static MatN Eye( const size_t n );
   void Eye(void);
   std::vector<double> GetVector(void) const { return m_mat; }
   void SetVector(const std::vector<double>& m) { m_mat = m; }

   virtual double operator() ( const size_t i );
   virtual double operator() ( const size_t row, const size_t col );
   virtual bool IsUnit( void ) const;

   virtual std::vector<double> GetMatrix( void ) const { return m_mat; }
   virtual void resize(const size_t n) {
      m_mat.resize(n);
      m_dim = n;
      m_rowdim = (size_t)(sqrt(double(n)));
   }

protected:
   size_t m_dim;
   std::vector<double> m_mat;
   size_t m_rowdim;


};

#endif

