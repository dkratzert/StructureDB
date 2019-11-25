#ifndef CELLINPUTDATA_H
#define CELLINPUTDATA_H

#include "LRL_Cell.h"

#include <algorithm>
#include <string>
#include <utility>

class CellInputData {
public:
   CellInputData(void);

   bool operator== (const std::string& s) const;
   std::string GetInputDataType(void) const { return m_inputDataType; }
   std::string GetStrCell(void) const { return m_strCell; }
   std::string GetLattice(void) const { return m_lattice; }
   LRL_Cell GetCell(void)           const { return m_cell; }
   std::pair<int, int> GetVarietyRange(void)           const { return m_varietyRange; }
   void SetInputDataType(const std::string& inputDataType) { m_inputDataType = inputDataType; }
   void SetLattice(const std::string& lattice) { m_lattice = lattice; }
   void SetStrCell(const std::string& s) { m_strCell = s; }
   void SetCell(const LRL_Cell& c) { m_cell = c; }
   void SetVarietyRange(const std::pair<int, int>& p) { m_varietyRange = p; }

protected:
   std::string m_input;
   std::string m_lattice;
   std::string m_inputDataType;
   std::string m_strCell;
   LRL_Cell        m_cell;
   std::pair<int, int> m_varietyRange;
};

#endif
