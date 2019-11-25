#ifndef SELLING_H
#define SELLING_H

class MatS6;
#include "S6.h"

class Selling {
public:
   template<typename TVEC>
   static bool Reduce(const TVEC& vin, TVEC& vout) {
      S6 s6out;
      const bool b = Selling::Reduce(TVEC(vin), s6out);
      vout = s6out;
      return b;
   }

   template<typename TVEC, typename TMAT>
   static bool Reduce(const TVEC& vin, TMAT& mReduce, TVEC& vout) {
      MatS6 m6Reduce;
      S6 s6out;
      const bool b = Selling::Reduce(S6(vin), m6Reduce, s6out);
      mReduce = m6Reduce;
      vout = s6out;
      return b;
   }

   static bool Reduce(const S6& in, S6& out);
   static bool Reduce(const S6& in, S6& out, const bool sellingFirst);

   static double MaxScalar(const S6& s6);

   static bool IsReduced(const S6& v);
   static bool IsReduced(const S6& v, const double delta);

   static void SetDebug(const bool b) { m_debugInstrument = b; }

   static size_t GetCycles() { 
      return m_ReductionCycleCount;
   }

   static void SetListSteps(const bool b) { listSteps = b; }

   static std::string GetName() { return "Selling"; }

private:
   static bool m_debugInstrument;
   static size_t m_ReductionCycleCount;
   static bool listSteps;
private:

};
#endif // SELLING_H