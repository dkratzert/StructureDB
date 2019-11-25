#pragma warning( disable : 4189) //   local variable is initialized but not referenced



#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "LRL_Cell.h"
#include "G6.h"
#include "S6.h"
#include "Selling.h"
#include "LRL_StringTools.h"


G6::G6(void) {
   m_vec.resize(6);
   m_dim = 6;
   m_valid = true;
}

G6::G6(const G6& v)
: G6()
{
   m_vec = v.m_vec;
   m_valid = v.m_valid;
}

G6::G6(const double v[6])
   : G6()
{
   m_vec[0] = v[0];
   m_vec[1] = v[1];
   m_vec[2] = v[2];
   m_vec[3] = v[3];
   m_vec[4] = v[4];
   m_vec[5] = v[5];
}

G6::G6(const LRL_Cell& c) {
   static const double pi = 4.0*atan(1.0);
   static const double twopi = 2.0*pi;

   m_vec.resize(6);
   m_valid = c.GetValid();
   m_vec[0] = c[0] * c[0];
   m_vec[1] = c[1] * c[1];
   m_vec[2] = c[2] * c[2];
   m_vec[3] = 2.0 * c[1] * c[2] * cos(c[3]);
   m_vec[4] = 2.0 * c[0] * c[2] * cos(c[4]);
   m_vec[5] = 2.0 * c[0] * c[1] * cos(c[5]);
   for (size_t i = 3; i<6; ++i) if (std::fabs(m_vec[i]) < 1.0E-10) m_vec[i] = 0.0;

   bool b1 = c.GetValid();
   bool b2 = c[3] < pi;
   bool b3 = c[4] < pi;
   bool b4 = c[5] < pi;
   bool b5 = (c[3] + c[4] + c[5])< twopi;
   m_valid = b1 && b2 && b3 && b4 && b5;
   m_vec.SetValid(m_valid);
}

G6::G6(const S6& ds)
   : G6()
{
   const double& p = ds[0];
   const double& q = ds[1];
   const double& r = ds[2];
   const double& s = ds[3];
   const double& t = ds[4];
   const double& u = ds[5];

   double& g1 = m_vec[0];
   double& g2 = m_vec[1];
   double& g3 = m_vec[2];
   double& g4 = m_vec[3];
   double& g5 = m_vec[4];
   double& g6 = m_vec[5];

   g1 = -s - r - q;
   g2 = -t - r - p;
   g3 = -u - q - p;
   g4 = 2.0*p;
   g5 = 2.0*q;
   g6 = 2.0*r;
   m_valid = g1 > 0.001 && g2 > 0.001 && g3 > 0.001 &&
      std::fabs(g4*g4 / g2 / g3) <= 4.0 &&
      std::fabs(g5*g5 / g1 / g3) <= 4.0 &&
      std::fabs(g6*g6 / g1 / g2) <= 4.0 &&
      ds.GetValid();
   m_vec.SetValid(m_valid);
 
}

G6::G6(const VecN& v) {
   m_dim = 0;
   m_valid = false;
   if (v.size() == 6) {
      m_dim = 6;
      m_valid = true;
      m_vec = v;
   }
}

G6::G6(const std::string& s) {
   m_vec.resize(6);
   m_vec = LRL_StringTools::FromString(s);
   m_dim = (size_t)(m_vec.size());
   m_valid = true;
   if (m_dim != 6) throw "bad dimension in G6 from a string";
}

G6::G6(const std::vector<double>& v) {
   m_vec.resize(6);
   m_dim = (size_t)(v.size());
   if (m_dim != 6) throw "bad dimension in G6 from a vector";
   m_vec = v;
   m_valid = true;
}

double G6::DistanceBetween(const G6& v1, const G6& v2) {
   double sum = 0.0;
   for (size_t i = 0; i < 6; ++i)  sum += (v1[i] - v2[i])*(v1[i] - v2[i]);
   return sqrt(sum);
}

G6& G6::operator= (const G6& v) {
   m_vec = v.m_vec;
   m_dim = v.m_dim;
   m_valid = v.m_valid;
   return *this;
}

G6& G6::operator= (const std::string& s)
{
   (*this) = G6(s);
   return *this;
}


G6& G6::operator= (const S6& v)
{
   (*this) = G6(v);
   return *this;
}


G6& G6::operator= (const LRL_Cell& v)
{
   static const double pi = 4.0*atan(1.0);
   static const double twopi = 2.0*pi;

   (*this) = G6(v);
   m_valid = v.GetValid() && v[3] < pi && v[4] < pi && v[5] < pi && (v[3] + v[4] + v[5])< twopi;
   return *this;
}

G6& G6::operator/= (const double d) {
   (*this).m_vec /= d;
   return *this;
}

G6& G6::operator*= (const double d) {
   (*this).m_vec *= d;
   return *this;
}

G6 G6::operator* (const double d) const {
   G6 g6(*this);
   g6.m_vec *= d;
   return g6;
}

G6 G6::operator/ (const double d) const {
   G6 g6(*this);
   g6.m_vec /= d;
   return g6;
}

G6 G6::operator- (void) const {
   G6 g6;
   g6.m_vec = -g6.m_vec;
   return *this; // unary
}

G6& G6::operator+= (const G6& g6) {
   for (size_t i = 0; i < g6.size(); ++i)
      m_vec[i] += g6.m_vec[i];
   return *this;
}

G6& G6::operator-= (const G6& g6) {
   for (size_t i = 0; i < g6.size(); ++i)
      m_vec[i] -= g6.m_vec[i];
   return *this;
}

bool G6::operator== (const G6& g6) const {
   return m_vec == g6.m_vec;
}

bool G6::operator!= (const G6& g6) const {
   return !((*this) == g6);
}

G6 G6::operator+ (const G6& v) const {
   G6 d;
   if ((*this).m_valid && v.m_valid) {
      d.m_vec = m_vec + v.m_vec;
   }
   else {
      d.m_valid = false;
   }
   return d;
}

G6 G6::operator- (const G6& v) const {
   G6 d;
   if ((*this).m_valid && v.m_valid) {
      d.m_vec = m_vec - v.m_vec;
   }
   else {
      d.m_valid = false;
   }
   return d;
}

G6 G6::InverseG6Vector(void) const {
   G6 v(*this);

   double a = sqrt(v[0]);
   double b = sqrt(v[1]);
   double c = sqrt(v[2]);
   double cosalpha = 0.5*v[3] / v[1] / v[2];
   double cosbeta = 0.5*v[4] / v[0] / v[2];
   double cosgamma = 0.5*v[5] / v[0] / v[1];

   const double sinalpha = sqrt(1.0 - cosalpha*cosalpha);
   const double sinbeta = sqrt(1.0 - cosbeta *cosbeta);
   const double singamma = sqrt(1.0 - cosgamma*cosgamma);

   const double volume = a*b*c*sqrt(1.0 - cosalpha*cosalpha - cosbeta*cosbeta - cosgamma*cosgamma) -
      2.0*cosalpha*cosbeta*cosgamma;

   const double astar = b*c*sinalpha / volume;
   const double bstar = a*c*sinbeta / volume;
   const double cstar = a*b*singamma / volume;

   const double cosalphastar = (cosbeta *cosgamma - cosalpha) / std::fabs(sinbeta *singamma);
   const double cosbetastar = (cosalpha*cosgamma - cosbeta) / std::fabs(sinalpha*singamma);
   const double cosgammastar = (cosalpha*cosbeta - cosgamma) / std::fabs(sinalpha*sinbeta);

   G6 vout;
   vout[0] = astar*astar;
   vout[1] = bstar*bstar;
   vout[2] = cstar*cstar;
   vout[3] = 2.0*bstar*cstar*cosalphastar;
   vout[4] = 2.0*astar*cstar*cosbetastar;
   vout[5] = 2.0*astar*bstar*cosgammastar;

   return vout;
}

G6 G6::sort(const G6& g6) {
   G6 g(g6);
   bool again = true;
   while (again) {
      if (g[1] < g[0]) {
         again = true;
         std::swap(g[0], g[1]);
         std::swap(g[3], g[4]);
      }
      else if (g[2] < g[1]) {
         again = true;
         std::swap(g[1], g[2]);
         std::swap(g[4], g[5]);
      }
      else
         again = false;
   }
   return g;
}

std::vector<std::pair<std::string, std::string> > G6::ClassifyVector(const double delta) const {
   std::vector<std::pair<std::string, std::string> > vss;
   const std::pair<std::string, std::string> noop = std::make_pair(".", "");
   const double g1 = (*this)[0];
   const double g2 = (*this)[1];
   const double g3 = (*this)[2];
   const double g4 = (*this)[3];
   const double g5 = (*this)[4];
   const double g6 = (*this)[5];

   if (std::fabs(g1 - g2) <= delta) {
      const std::string s = "1";
      vss.push_back(std::make_pair("1", "1"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g2 - g3) <= delta) {
      const std::string s = "2";
      vss.push_back(std::make_pair("2", "2"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g4) <= delta) {
      const std::string s = "3";
      vss.push_back(std::make_pair("3", "3"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g5) <= delta) {
      const std::string s = "4";
      vss.push_back(std::make_pair("4", "4"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g6) <= delta) {
      const std::string s = "5";
      vss.push_back(std::make_pair("5", "5"));
   }
   else {
      vss.push_back(noop);
   }

   if ((std::fabs(g2 - g4) <= delta) && ((g5 - g6) >= -delta)) {
      const std::string s = "6";
      vss.push_back(std::make_pair("6", "6"));
   }
   else {
      vss.push_back(noop);
   }

   if ((std::fabs(g2 - g4) <= delta) + ((g5 - g6) <= delta)) {
      const std::string s = "7";
      vss.push_back(std::make_pair("7", "7"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g2 + g4) <= delta) {
      const std::string s = "8";
      vss.push_back(std::make_pair("8", "8"));
   }
   else {
      vss.push_back(noop);
   }

   if ((std::fabs(g1 - g5)) <= delta && ((g4 - g6) >= -delta)) {
      const std::string s = "9";
      vss.push_back(std::make_pair("9", "9"));
   }
   else {
      vss.push_back(noop);
   }

   if ((std::fabs(g1 - g5)) <= delta && ((g4 - g6) <= delta)) {
      const std::string s = "A";
      vss.push_back(std::make_pair("A", "A"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g1 + g5) <= delta) {
      const std::string s = "B";
      vss.push_back(std::make_pair("B", "B"));
   }
   else {
      vss.push_back(noop);
   }

   if ((std::fabs(g1 - g6) <= delta) && ((g4 - g5) >= -delta)) {
      const std::string s = "C";
      vss.push_back(std::make_pair("C", "C"));
   }
   else {
      vss.push_back(noop);
   }

   if ((std::fabs(g1 - g6) <= delta) && ((g4 - g5) <= delta)) {
      const std::string s = "D";
      vss.push_back(std::make_pair("D", "D"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g1 + g6) <= delta) {
      const std::string s = "E";
      vss.push_back(std::make_pair("E", "E"));
   }
   else {
      vss.push_back(noop);
   }

   if (std::fabs(g1 + g2 + g3 + g4 + g5 + g6 - g3) <= delta) {
      const std::string s = "F";
      vss.push_back(std::make_pair("F", "F"));
   }
   else {
      vss.push_back(noop);
   }

   return vss;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
std::pair<int, std::string> G6::IdentifyNearbyBoundaries(const G6& v, const double cutoff)
/*-------------------------------------------------------------------------------------*/
{
   std::string s1, s2;

   if (std::fabs(v[0] - v[1]) < cutoff)
   {
      s1 += '1';
      s2 += '1';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[1] - v[2]) < cutoff)
   {
      s1 += '2';
      s2 += '2';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[3]) < cutoff)
   {
      s1 += '3';
      s2 += '3';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[4]) < cutoff)
   {
      s1 += '4';
      s2 += '4';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[5]) < cutoff)
   {
      s1 += '5';
      s2 += '5';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[1] - v[3]) < cutoff &&
      (v[5] - v[4]) < cutoff)
   {
      s1 += '6';
      s2 += '6';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[1] - v[3]) < cutoff &&
      (v[4] - v[5]) < cutoff)
   {
      s1 += '7';
      s2 += '7';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[1] + v[3]) < cutoff)
   {
      s1 += '8';
      s2 += '8';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[0] - v[4]) < cutoff &&
      (v[5] - v[3]) < cutoff)
   {
      s1 += '9';
      s2 += '9';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[0] - v[4]) < cutoff &&
      (v[3] - v[5]) < cutoff)
   {
      s1 += 'A';
      s2 += 'A';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[0] + v[4]) < cutoff)
   {
      s1 += 'B';
      s2 += 'B';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[0] - v[5]) < cutoff &&
      (v[4] - v[3]) < cutoff)
   {
      s1 += 'C';
      s2 += 'C';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[0] - v[5]) < cutoff &&
      (v[3] - v[4]) < cutoff)
   {
      s1 += 'D';
      s2 += 'D';
   }
   else
   {
      s2 += '.';
   }

   if (std::fabs(v[0] + v[5]) < cutoff)
   {
      s1 += 'E';
      s2 += 'E';
   }
   else
   {
      s2 += '.';
   }

   if (fabs(v[0] + v[1] + v[3] + v[4] + v[5]) < cutoff)
   {
      s1 += 'F';
      s2 += 'F';
   }
   else
   {
      s2 += '.';
   }

   return std::make_pair(-1, (s1.empty()) ? "" : "(" + s1 + ") \t(" + s2 + ")");
}

static bool TestCase1(const G6& g) { return g[0] == g[1]; }
static bool TestCase2(const G6& g) { return g[1] == g[2]; }
static bool TestCase3(const G6& g) { return g[3] == 0; }
static bool TestCase4(const G6& g) { return g[4] == 0; }
static bool TestCase5(const G6& g) { return g[5] == 0; }
static bool TestCase6(const G6& g) { return (g[1] == g[3]) && (g[4] >= g[5]); }
static bool TestCase7(const G6& g) { return (g[1] == g[3]) && (g[4] < g[5]); }
static bool TestCase8(const G6& g) { return g[1] == -g[3]; }
static bool TestCase9(const G6& g) { return (g[0] == g[4]) && (g[3] >= g[5]); }
static bool TestCaseA(const G6& g) { return (g[0] == g[4]) && (g[3] < g[5]); }
static bool TestCaseB(const G6& g) { return g[0] == -g[4]; }
static bool TestCaseC(const G6& g) { return (g[0] == g[5]) && (g[3] >= g[4]); }
static bool TestCaseD(const G6& g) { return (g[0] == g[5]) && (g[3] < g[4]); }
static bool TestCaseE(const G6& g) { return g[0] == -g[5]; }
static bool TestCaseF(const G6& g) {
   return (g[0] + g[1] + g[2] + g[3] + g[4] + g[5]) == g[2];
}

bool G6::TestBoundaryDetection(void) {
   bool b1 = true;
   b1 = b1 && TestCase1(G6("1 1 0 0 0 0"));
   b1 = b1 && TestCase2(G6("0 1 1 0 0 0"));
   b1 = b1 && TestCase3(G6("1 2 3 0 5 6"));
   b1 = b1 && TestCase4(G6("1 2 3 4 0 6"));
   b1 = b1 && TestCase5(G6("1 2 3 4 5 0"));
   b1 = b1 && TestCase6(G6("0 1 0 1 9 3"));
   b1 = b1 && TestCase7(G6("0 1 0 1 3 9"));
   b1 = b1 && TestCase8(G6("0 1 0 -1 0 0"));
   b1 = b1 && TestCase9(G6("1 0 0 3 1 2"));
   b1 = b1 && TestCaseA(G6("1 0 0 2 1 3"));
   b1 = b1 && TestCaseB(G6("1 0 0 0 -1 0"));
   b1 = b1 && TestCaseC(G6("1 0 0 9 2 1"));
   b1 = b1 && TestCaseD(G6("1 0 0 2 9 1"));
   b1 = b1 && TestCaseE(G6("1 0 0 0 0 -1"));
   b1 = b1 && TestCaseF(G6("1 3 0 5 7 -16"));

   std::vector<std::pair<std::string, std::string> > s1 = G6("1 1 0 0 0 0"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s2 = G6("0 1 1 0 0 0"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s3 = G6("1 2 3 0 5 6"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s4 = G6("1 2 3 4 0 6"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s5 = G6("1 2 3 4 5 0"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s6 = G6("0 1 0 1 9 3"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s7 = G6("0 1 0 1 3 9"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s8 = G6("0 1 0 -1 0 0" ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > s9 = G6("1 0 0 3 1 2"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > sA = G6("1 0 0 2 1 3"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > sB = G6("1 0 0 0 -1 0" ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > sC = G6("1 0 0 9 2 1"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > sD = G6("1 0 0 2 9 1"  ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > sE = G6("1 0 0 0 0 -1" ).ClassifyVector(0.0);
   std::vector<std::pair<std::string, std::string> > sF = G6("1 3 0 5 7 -16").ClassifyVector(0.0);

   return b1;
}

std::ostream& operator<< (std::ostream& o, const G6& v) {
   std::streamsize oldPrecision = o.precision();
   o << std::fixed << std::setprecision(3);
   for (size_t i = 0; i < v.size(); ++i)
      o << std::setw(9) << v[i] << " ";
   o << std::setprecision(oldPrecision);
   o.unsetf(std::ios::floatfield);
   return o;
}

G6 operator*(const double d, const G6& v) {
   return v*d;
}
