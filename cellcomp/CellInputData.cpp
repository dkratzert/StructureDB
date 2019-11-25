


#include "CellInputData.h"

CellInputData::CellInputData( void ) 
: m_lattice( "")
, m_inputDataType("")
, m_strCell("")
{ }

bool CellInputData::operator== (const std::string& s) const {
   return m_lattice == s;
}
