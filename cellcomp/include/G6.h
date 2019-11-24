#ifndef G6_H
#define G6_H

#include "VecN.h"

class LRL_Cell;
class S6;

#include <ostream>
#include <string>

class G6 {
public:
   friend std::ostream& operator<< (std::ostream&, const G6&);
   friend G6 operator* (const double d, const G6& v);
   friend class S6;

   G6(void);
   G6(const G6& v);
   G6(const double v[6]);
   G6(const LRL_Cell& c);
   G6(const VecN& v);
   G6(const S6& ds);
   G6(const std::string& s);
   G6(const std::vector<double>& v);
   ~G6(void) {}

   G6& operator= (const G6& v);
   G6& operator= (const std::string& s);
   G6& operator= (const S6& v);
   G6& operator= (const LRL_Cell& v);
   G6& operator/= (const double d);
   G6& operator*= (const double d);
   G6& operator+= (const G6& v);
   G6& operator-= (const G6& v);
   G6 operator+ (const G6& v) const;
   G6 operator- (const G6& v) const;
   G6 operator* (const double d) const;
   G6 operator/ (const double d) const;
   G6 operator- (void) const; // unary
   bool operator== (const G6& g6) const;
   bool operator!= (const G6& g6) const;


   void SetVector( const std::vector<double>& v ) { m_vec = v; }
   double operator[](const size_t n) const { return m_vec[n]; }
   double& operator[](const size_t n) { return m_vec[n]; }
   double DistanceBetween(const G6& v1, const G6& v2);
   std::vector<std::pair<std::string, std::string> > ClassifyVector(const double delta) const;
   bool TestBoundaryDetection(void);
   size_t size(void) const { return 6; }
   double norm() const { return m_vec.norm(); }
   double norm(const G6& g6) const { return g6.norm(); }
   double Norm() const { return m_vec.norm(); }
   double Norm(const G6& g6) const { return g6.norm(); }

   G6 InverseG6Vector(void) const;
   G6 sort(const G6& g6);
   std::vector<double> GetVector(void) const { return m_vec.GetVector(); }
   double* data() const { return const_cast<double*>(m_vec.data()); }
   bool GetValid(void) const { return m_valid; }
   bool IsValid(void) const { return m_valid; }
   void SetValid(const bool b) { m_valid = b; }

   double at(const size_t n) const { return m_vec[n]; }

   static std::string GetName(void) { return "G6"; }

   static std::pair<int, std::string> IdentifyNearbyBoundaries(const G6& v, const double cutoff);

private:
   VecN m_vec;
   size_t m_dim;
   bool m_valid;
};

#endif

/*
Niggli reduction consists of converting the original cell to a primitive one and then
alternately applying two operations: conversion to standard presentation and reduction [?]
[?, Table 1]. The convention for meeting the combined Buerger and Niggli conditions is
based on increasingly restrictive layers of constraints:
If g1 < g2 < g3, |g4| < g2, |g5| < g1, |g6| < g1 and either g{4,5,6} > 0 or g{4,5,6}<=0 then
we have a Niggli-reduced cell, and we are done.
The remaining conditions are imposed when any of the above inequalities becomes
an equality or the elements of gf4;5;6g are not consistently all strictly positive or are not
consistently all less than or equal to zero.
The full set of combined Buerger and Niggli conditions in addition to those for the cell
edge lengths being minimal is:

require 0<=g1<=g2<=g3
if g2 = g3; then require |g5|<=|g6|
require g4 > 0 and g5 > 0 and g6 > 0
or require g4<=0 and g5<=0 and g6<=0
require |g4|<=g2
require |g5|<=g1
require |g6|<=g1
require |g3|<=g1 + g2 + g3 + g4 + g5 + g6
if g4 = g2; then require g6<=2g5
if g5 = g1; then require g6<=2g4
if g6 = g1; then require g5<=2g4
if g4 = g2; then require g6 = 0
if g5 = g1; then require g6 = 0
if g6 = g1; then require g5 = 0
if g3 = g1 + g2 + g3 + g4 + g5 + g6; then require 2g1 + 2g5 + g6<=0
The G6 transformations associated with each of these steps are enumerated in [?]. Appli-
cation of these operations must be repeated until all are satisfed.
*/