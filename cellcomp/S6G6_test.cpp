// S6G6.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include "LatticeConverter.h"
#include "LRL_Cell.h"
#include "LRL_ReadLatticeData.h"
#include "NCDist.h"
#include "S6.h"
#include "CS6Dist.h"


std::string Letters( void ) {
   return "V,G,D,S,P,A,B,C,I,F,R,C3,G6,S6,B4,D7,H";
}

std::vector<LRL_ReadLatticeData> GetInputCells( void ) {
   const std::string letters = Letters( );
   std::string lattice;
   std::vector<LRL_ReadLatticeData> cellDataList;
   LRL_ReadLatticeData rcd;
   while (lattice != "EOF") {
      rcd.read( );
      lattice = rcd.GetLattice( );
      if ((!lattice.empty( )) && (letters.find( static_cast<char>(toupper( lattice[0] )) ) != std::string::npos))
         cellDataList.push_back( rcd );
   }

   return cellDataList;
}

std::vector<S6> GetInputSellingReducedVectors( const std::vector<LRL_ReadLatticeData>& input ) {
   std::vector<S6> v;
   LatticeConverter converter;


   for (size_t i = 0; i < input.size( ); ++i) {
      const S6 s6 = converter.SellingReduceCell( input[i].GetLattice( ), input[i].GetCell( ) );
      v.push_back( s6 );
   }
   return v;
}


int main()
{
   const std::vector<LRL_ReadLatticeData> input = GetInputCells( );

   const std::vector<S6> vLat = GetInputSellingReducedVectors( input );

   for (size_t lat = 0; lat < vLat.size( ); ++lat) {
      if (vLat[lat].IsValid( )) {
         const LRL_Cell cell1( vLat[lat] );
         const S6 a( cell1 );
         for (size_t i = lat + 1; i < vLat.size( ); ++i) {
            if (vLat[i].IsValid( )) {
               std::cout << std::endl;
               const LRL_Cell cell2( vLat[i] );
               const S6 b( cell2 );
               std::cout << lat << "  " << cell1 << std::endl;
               std::cout << i << "  " << cell2 << std::endl;
               std::cout << lat << "  " << a << std::endl;
               std::cout << i << "  " << b << std::endl;
               std::cout << lat << "  " << i << "  " << CS6Dist( a.data( ), b.data( ) ) << std::endl;
               std::cout << "Angstrom distance  " << 0.1*sqrt( CS6Dist( a.data( ), b.data( ) ) );
            }
         }
      }
   }

   const int  i19191 = 19191;
}
