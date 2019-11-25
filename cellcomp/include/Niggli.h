#ifndef REDUCER_H
#define REDUCER_H

#include "G6.h"


class Niggli {
public:
   ~Niggli(void);
   static bool Reduce(const G6& vi, G6& vout);
   static bool IsNiggli(const G6& s);
   static void MKnormWithoutMatrices(const G6& vi, G6& vout, const double delta);
   static bool ReduceWithoutMatrices(const G6& vi, G6& vout, const double delta);

   static std::string GetName() { return "Niggli"; }
   static size_t GetCycles() { return m_ReductionCycleCount; }

private:
   static size_t m_ReductionCycleCount;

   // at least for now, all functions are static, and there is no member data
   // forbid constructor
   Niggli(void);

public:

};

#endif // REDUCER_H

