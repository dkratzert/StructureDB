#pragma warning( disable : 4100) //  unreferenced formal parameter




#include "Selling.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <utility>


bool Selling::m_debugInstrument = false;
size_t Selling::m_ReductionCycleCount = 0;

std::vector< S6(*)(const S6&)> FillReduceFunctionArray() {
   static std::vector< S6(*)(const S6&)> vf;
   if (vf.empty()) {
      // THESE ARE TRANSFORMATIONS IN S6 (I THINK), NOT NOT NOT IN G6 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      vf.push_back(&S6::Reduce11);  // g or p
      vf.push_back(&S6::Reduce21);  // h or q
      vf.push_back(&S6::Reduce31);  // k or r
      vf.push_back(&S6::Reduce41);  // l or s
      vf.push_back(&S6::Reduce51);  // m or t
      vf.push_back(&S6::Reduce61);  // n or u
   }
   return vf;
}

bool Selling::Reduce(const S6& in, S6& out, const bool sellingFirst) {
   // dummy function to have the same signature as one in Niggli
   return Reduce(in, out);
}

bool Selling::Reduce(const S6& in, S6& out) {
   static const std::vector< S6(*)(const S6&)> reductionFunctions = FillReduceFunctionArray();
   size_t maxIndex = INT_MAX;
   m_ReductionCycleCount = 0;
   out = in;
   if (out[0] <= 0.0 && out[1] <= 0.0 && out[2] <= 0.0 && out[3] <= 0.0 && out[4] <= 0.0 && out[5] <= 0.0) return true;

   while (S6::CountPositive(out) != 0) {
      double maxScalar = -DBL_MAX;

      for (size_t i = 0; i < 6; ++i) {
         if (out[i] > maxScalar) {
            maxIndex = i;
            maxScalar = out[i];
         }
      }
      
      out = reductionFunctions[maxIndex](out);
      ++m_ReductionCycleCount;
      if (m_ReductionCycleCount > 1000 || S6::NegativeSumOfScalars(out) < 0.0 ) return false;
   }
   return true;
}

   double Selling::MaxScalar(const S6& s6) {
      const std::vector<double>& v(s6.GetVector());
      return *std::max_element(v.begin(), v.end());
   }

   bool Selling::IsReduced(const S6& s6) {
      return (MaxScalar(s6) <= 0.0);
   }

   bool Selling::IsReduced(const S6& s6, const double delta) {
      return (MaxScalar(s6) <= delta);
   }
