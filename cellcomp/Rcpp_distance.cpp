#include <R.h>
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace RcppParallel;


#include <vector>

#include "CS6Dist.h"
#include "G6.h"
#include "LatticeConverter.h"
#include "LRL_Cell.h"
#include "NCDist.h"
#include "S6.h"

extern "C" SEXP rcpp_distance( SEXP choice, SEXP lat1_, SEXP a1_, SEXP b1_, SEXP c1_,
   SEXP alpha1_, SEXP beta1_, SEXP gamma1_,
   SEXP lat2_, SEXP a2_, SEXP b2_, SEXP c2_,
   SEXP alpha2_, SEXP beta2_, SEXP gamma2_ ) {
   G6 prim1, prim2;

   const char typeDistance = toupper(Rcpp::as<std::string>( choice )[0]);

   double dcell1[6],dcell2[6];
   LRL_Cell cell1, cell2;

   std::string lat1 = Rcpp::as<std::string>( lat1_ );
   dcell1[0] = Rcpp::as<double>( a1_ );
   dcell1[1] = Rcpp::as<double>( b1_ );
   dcell1[2] = Rcpp::as<double>( c1_ );
   dcell1[3] = Rcpp::as<double>( alpha1_ );
   dcell1[4] = Rcpp::as<double>( beta1_ );
   dcell1[5] = Rcpp::as<double>( gamma1_ );
   cell1 = LRL_Cell(dcell1[0],dcell1[1],dcell1[2],dcell1[3],dcell1[4],dcell1[5]);

   std::string lat2 = Rcpp::as<std::string>( lat2_ );
   dcell2[0] = Rcpp::as<double>( a2_ );
   dcell2[1] = Rcpp::as<double>( b2_ );
   dcell2[2] = Rcpp::as<double>( c2_ );
   dcell2[3] = Rcpp::as<double>( alpha2_ );
   dcell2[4] = Rcpp::as<double>( beta2_ );
   dcell2[5] = Rcpp::as<double>( gamma2_ );
   cell2 = LRL_Cell(dcell2[0],dcell2[1],dcell2[2],dcell2[3],dcell2[4],dcell2[5]);

   double answer = -1.0;

   if (typeDistance == 'S') {
      const S6 prim1 = LatticeConverter::SellingReduceCell( lat1, cell1 );
      const S6 prim2 = LatticeConverter::SellingReduceCell( lat2, cell2 );

      answer = 0.1*std::sqrt( CS6Dist( prim1.data( ), prim2.data( ) ) );
    /*Rcpp::Rcout<< "cell1: ["<<cell1[0]<<","<<cell1[1]<<","<<cell1[2]<<"," 
      <<cell1[3]<< ","<<cell1[4]<<","<<cell1[5]<<"]"<<std::endl;
      Rcpp::Rcout<< "s6primred1: ["<<prim1[0]<<","<<prim1[1]<<","<<prim1[2]<<"," 
      <<prim1[3]<< ","<<prim1[4]<<","<<prim1[5]<<"]"<<std::endl;
      Rcpp::Rcout<< "cell2: ["<<cell2[0]<<","<<cell2[1]<<","<<cell2[2]<<"," 
      <<cell2[3]<< ","<<cell2[4]<<","<<cell2[5]<<"]"<<std::endl;
      Rcpp::Rcout<< "s6primred2: ["<<prim2[0]<<","<<prim2[1]<<","<<prim2[2]<<"," 
      <<prim2[3]<< ","<<prim2[4]<<","<<prim2[5]<<"]"<<std::endl;
      Rcpp::Rcout<<"answer: "<<answer<<std::endl; */
   } else if (typeDistance == 'G' || typeDistance == 'N') {
      const G6 prim1 = LatticeConverter::NiggliReduceCell( lat1, cell1 );
      const G6 prim2 = LatticeConverter::NiggliReduceCell( lat2, cell2 );

      answer = 0.1*std::sqrt( NCDist( prim1.data( ), prim2 .data( ) ) );
    /*Rcpp::Rcout<< "cell1: ["<<cell1[0]<<","<<cell1[1]<<","<<cell1[2]<<"," 
      <<cell1[3]<< ","<<cell1[4]<<","<<cell1[5]<<"]"<<std::endl;
      Rcpp::Rcout<< "g6primred1: ["<<prim1[0]<<","<<prim1[1]<<","<<prim1[2]<<"," 
      <<prim1[3]<< ","<<prim1[4]<<","<<prim1[5]<<"]"<<std::endl;
      Rcpp::Rcout<< "cell2: ["<<cell2[0]<<","<<cell2[1]<<","<<cell2[2]<<"," 
      <<cell2[3]<< ","<<cell2[4]<<","<<cell2[5]<<"]"<<std::endl;
      Rcpp::Rcout<< "g6primred2: ["<<prim2[0]<<","<<prim2[1]<<","<<prim2[2]<<"," 
      <<prim2[3]<< ","<<prim2[4]<<","<<prim2[5]<<"]"<<std::endl;
      Rcpp::Rcout<<"answer: "<<answer<<std::endl;*/
   }

   return Rcpp::wrap( answer );
}
