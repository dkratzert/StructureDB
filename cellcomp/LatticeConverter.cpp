#include "G6.h"
#include "LRL_Cell.h"
#include "MatG6.h"
#include "S6.h"
#include "LatticeConverter.h"
#include "LRL_ToString.h"
#include "Niggli.h"
#include "Selling.h"

#include <string>


void LatticeConverter::TextOutput(const std::string& label, const std::string& lattice, const LRL_Cell& cell) const {
   std::cout << label << std::endl;
   std::cout << "lattice " << lattice << std::endl;
   std::cout << "LRL_Cell  " << LRL_ToString(cell) << std::endl;
   std::cout << "G6 " << LRL_ToString(G6(cell)) << std::endl;
   std::cout << "Scalars " << LRL_ToString(S6(cell)) << std::endl;
}

LatticeConverter::LatticeConverter(const eOutputType type)
   : m_OutputType(type)
{

}

LRL_Cell LatticeConverter::MakePrimitiveCell(const std::string& lattice, const LRL_Cell& cell) {
   const G6 g6 = G6(cell);
   const MatG6 mLattice = LRL_Cell::LatSymMatG6(lattice, cell);
   G6 redVec;
   return LRL_Cell(mLattice*g6);
}

void LatticeConverter::SetOutputMaxima(void) { m_OutputType = emaxima; }
void LatticeConverter::SetOutputText  (void) { m_OutputType = etext; }

LRL_Cell LatticeConverter::NiggliReduceCell(const std::string& lattice, const LRL_Cell& cell) {
   const G6 g6 = G6(cell);
   const MatG6 mLattice = LRL_Cell::LatSymMatG6(lattice, cell);
   MatG6 m66;
   G6 redVec;
   const bool b = Niggli::Reduce(mLattice*g6, redVec);
   if (b) {
      return LRL_Cell(redVec);;
   }
   else {
      return LRL_Cell();
   }
}

LRL_Cell LatticeConverter::SellingReduceCell(const std::string& lattice, const LRL_Cell& cell) {
   const G6 g6 = G6(cell);
   const S6 tempS6 = g6;
   const MatG6 mLattice = LRL_Cell::LatSymMatG6(lattice, cell);
   S6 redVec;
   const bool b = Selling::Reduce(S6(mLattice*g6), redVec);
   if (b) {
      return LRL_Cell(redVec);;
   }
   else {
      return LRL_Cell();
   }
}
