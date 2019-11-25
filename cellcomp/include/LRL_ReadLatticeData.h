#ifndef LRL_READLATTICEDATA_H
#define LRL_READLATTICEDATA_H

#include "LRL_Cell.h"
#include "CellInputData.h"

class LRL_ReadLatticeData : public CellInputData {
public:
   LRL_ReadLatticeData read();
   void CellReader(const std::string& s);
   void CellReader(const std::string& lattice, const std::string& cell);
   bool IsLatticeName(const std::string inputName, const std::vector<std::string>& nameList);
   bool IsLatticeName(const std::vector<std::string>& nameList, const std::string inputName);

   static std::vector<CellInputData> ReadAllLatticeData();
   static std::vector<CellInputData> ReadAllLatticeDataAndMakePrimitive();
};

#endif // LRL_READLATTICEDATA_H
