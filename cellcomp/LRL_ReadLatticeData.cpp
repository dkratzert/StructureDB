


#include "LRL_Cell.h"
#include "G6.h"
#include "S6.h"
#include "LatticeConverter.h"
#include "LRL_ReadLatticeData.h"
#include "LRL_StringTools.h"
#include "LRL_ToString.h"
#include "LRL_StringTools.h"


#include <algorithm>
#include <complex>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>

void LRL_ReadLatticeData::CellReader(const std::string& lattice, const std::string& cell) {
   CellReader(lattice + " " + cell);
}

static int folseed = 19192;

//static const std::vector<size_t> vAvoid(avoidList, avoidList + sizeof(avoidList) / sizeof(avoidList[0]));
std::string G6_names[] = { "G ", "V ", "G6 " };
std::vector<std::string> vG6_names(G6_names, G6_names + sizeof(G6_names) / sizeof(G6_names[0]));

std::string S6_names[] = { "S ", "S6 " };
std::vector<std::string> vS6_names(S6_names, S6_names + sizeof(S6_names) / sizeof(S6_names[0]));

std::string lattice_names[] = { "A ", "B ", "C ", "P ", "R ", "F ", "I ", "H " };
std::vector<std::string> vlattice_names(lattice_names, lattice_names + sizeof(lattice_names) / sizeof(lattice_names[0]));

bool LRL_ReadLatticeData::IsLatticeName(const std::string inputName, const std::vector<std::string>& nameList) {
   return std::find(nameList.begin(), nameList.end(), inputName) != nameList.end();
}

bool LRL_ReadLatticeData::IsLatticeName(const std::vector<std::string>& nameList, const std::string inputName) {
   return std::find(nameList.begin(), nameList.end(), inputName) != nameList.end();
}

void LRL_ReadLatticeData::CellReader(const std::string& s) {
   bool valid = false;
   m_input = s;
   if ((LRL_StringTools::strToupper(s.substr(0, 3)) == std::string("END"))) {
      m_lattice = "EOF";
   }
   else {
      std::istringstream iss(s);
      iss >> m_inputDataType;
      m_inputDataType = LRL_StringTools::strToupper(m_inputDataType + " ");
     if ( IsLatticeName(m_inputDataType, vG6_names) ) {
         G6 v6;
         for (size_t i = 0; i < 6; ++i)
            iss >> v6[i];
         m_cell = v6;
         valid = m_cell.GetValid();
      }
      else if (IsLatticeName(m_inputDataType, vS6_names)) {
         S6 e;
         for (size_t i = 0; i < 6; ++i)
            iss >> e[i];
         valid = true;
         m_cell = e;
      }
      else if (IsLatticeName(m_inputDataType, vlattice_names)) {
         for (size_t i = 0; i < 6; ++i)
            iss >> m_cell[i];
         for (size_t i = 3; i < 6; ++i)
            m_cell[i] *= 4.0 * atan(1.0) / 180.0;
         valid = true;
         m_cell.SetValid(true);
      }
      else {
      }

      if (valid && m_cell.GetValid()) {
         iss >> m_lattice;
         m_lattice = m_inputDataType[0];
         m_cell.SetValid(true);
         if (LRL_StringTools::strToupper(m_lattice) == "END" || LRL_StringTools::strToupper(m_lattice) == "EOF") {
            m_lattice = "EOF";
            return;
         }
      }
      else {
         std::cout << "input line rejected, invalid" << std::endl;
         m_lattice = "";
      }
   }
}

std::vector<CellInputData> LRL_ReadLatticeData::ReadAllLatticeData() {
   std::vector<CellInputData> celldata;
   std::string lattice("");
   LRL_ReadLatticeData rcd;
   while (lattice != "EOF") {
      rcd.read();
      lattice = "EOF";
      lattice = rcd.GetLattice();
      if (lattice != "EOF" && !lattice.empty()) {
         rcd.SetVarietyRange(std::make_pair(0, 23));
         celldata.push_back(rcd);
      }
   }
   return celldata;
}

std::vector<CellInputData> LRL_ReadLatticeData::ReadAllLatticeDataAndMakePrimitive(){
   std::vector<CellInputData> cellData = ReadAllLatticeData();
   for ( size_t i=0; i<cellData.size(); ++i ) {
      if ( cellData[i].GetLattice() != "P") {
         LRL_Cell cell = LatticeConverter::MakePrimitiveCell(cellData[i].GetLattice(), cellData[i].GetCell());
         cellData[i].SetCell(cell);
         cellData[i].SetLattice("P");
      }
   }
   return cellData;
}

LRL_ReadLatticeData LRL_ReadLatticeData::read(void) {
   std::getline(std::cin, m_strCell);
   if (std::cin && (LRL_StringTools::strToupper(m_strCell.substr(0, 3)) != std::string("END"))) {
      CellReader(m_strCell);
      return *this;
   }
   else {
      m_lattice = "EOF";
      return *this;
   }
}
