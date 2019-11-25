//
//  CS6Dist.c
//
//
//  Created by Herbert J. Bernstein and Lawrence C. Andrews on 3/04/17.
//
//

#ifndef CS6DIST_H
#define CS6DIST_H


/* #define S6DIST_DEBUG */
#define S6DIST_NO_OUTER_PASS  

#include <math.h>
#include <float.h>

#ifdef S6DIST_DEBUG
static int CS6M_changed = 0;
static double oldvalue;
#include <stdio.h>
#define CSM_report_double(prolog,value,epilog) \
oldvalue=value; fprintf(stderr,"%s%g%s",prolog,value,epilog);
#define CSM_report_integer(prolog,value,epilog) \
fprintf(stderr,"%s%d%s",prolog,value,epilog);
#define CSM_report_double_if_changed(prolog,value,epilog) \
CS6M_changed=0; if (fabs(value-oldvalue)>1.e-8*(fabs(value)+fabs(oldvalue)+1.e-12)) {oldvalue=value; CS6M_changed=1; fprintf(stderr,"%s%g%s",prolog,value,epilog);}
#define CS6M_also_if_changed_report(prolog,value,epilog) \
if(CS6M_changed) {fprintf(stderr,"%s%s%s",prolog,value,epilog);}
#define CS6M_also_if_changed_report_integer(prolog,value,epilog) \
if(CS6M_changed) {fprintf(stderr,"%s%d%s",prolog,value,epilog);}
#define CS6M_also_if_changed_report_double(prolog,value,epilog) \
if(CS6M_changed) {fprintf(stderr,"%s%g%s",prolog,value,epilog);}
#define CS6M_also_if_changed_report_double_vector(prolog,value,epilog) \
if(CS6M_changed) {fprintf(stderr,"%s[%g, %g, %g, %g, %g, %g]%s",prolog,value[0],value[1],value[2],value[3],value[4],value[5],epilog);}
#define CSM_report_double_vector(prolog,value,epilog) \
{fprintf(stderr,"%s[%g, %g, %g, %g, %g, %g]%s",prolog,value[0],value[1],value[2],value[3],value[4],value[5],epilog);}
#else
#define CSM_report_double(prolog,value,epilog)
#define CSM_report_integer(prolog,value,epilog)
#define CSM_report_double_if_changed(prolog,value,epilog)
#define CS6M_also_if_changed_report(prolog,value,epilog)
#define CS6M_also_if_changed_report_integer(prolog,value,epilog)
#define CS6M_also_if_changed_report_double(prolog,value,epilog)
#define CS6M_also_if_changed_report_double_vector(prolog,value,epilog)
#define CSM_report_double_vector(prolog,value,epilog)
#endif


#define CS6M_min(a,b) (((a)<(b))?(a):(b))
#define CS6M_max(a,b) (((a)<(b))?(b):(a))

#define NS6BND 6

#define S6P_1 0
#define S6P_2 1
#define S6P_3 2
#define S6P_4 3
#define S6P_5 4
#define S6P_6 5

#define NS6BDPRJ 6


/* S6Dist Boundary projectors */

static double s6prj_perp[NS6BDPRJ][36] = {
   /*prj[S6P_1]  ,1,*/
   { 1.,0.,0.,0.,0.,0.,
   0.,0.,0.,0.,0.,0.,
   0.,0.,0.,0.,0.,0.,
   0.,0.,0.,0.,0.,0.,
   0.,0.,0.,0.,0.,0.,
   0.,0.,0.,0.,0.,0. },

   /*prj[S6P_2]  ,2,*/
{ 0.,0.,0.,0.,0.,0.,
0.,1.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0. },

/*prj[S6P_3]  ,3,*/
{ 0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,1.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0. },

/*prj[S6P_4]  ,4,*/
{ 0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,1.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0. },

/*prj[S6P_5]  ,5,*/
{ 0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,1.,0.,
0.,0.,0.,0.,0.,0. },

/*prj[S6P_6]  ,6,*/
{ 0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,1. } };


static double s6prj[NS6BDPRJ][36] = {
   /*prj_perp[S6P_1]  ,1,*/
   { 0.,0.,0.,0.,0.,0.,
   0.,1.,0.,0.,0.,0.,
   0.,0.,1.,0.,0.,0.,
   0.,0.,0.,1.,0.,0.,
   0.,0.,0.,0.,1.,0.,
   0.,0.,0.,0.,0.,1. },

   /*prj_perp[S6P_2]  ,2,*/
{ 1.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,1.,0.,0.,0.,
0.,0.,0.,1.,0.,0.,
0.,0.,0.,0.,1.,0.,
0.,0.,0.,0.,0.,1. },

/*prj_perp[S6P_3]  ,3,*/
{ 1.,0.,0.,0.,0.,0.,
0.,1.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,1.,0.,0.,
0.,0.,0.,0.,1.,0.,
0.,0.,0.,0.,0.,1. },

/*prj_perp[S6P_4]  ,4,*/
{ 1.,0.,0.,0.,0.,0.,
0.,1.,0.,0.,0.,0.,
0.,0.,1.,0.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,1.,0.,
0.,0.,0.,0.,0.,1. },

/*prj_perp[S6P_5]  ,5,*/
{ 1.,0.,0.,0.,0.,0.,
0.,1.,0.,0.,0.,0.,
0.,0.,1.,0.,0.,0.,
0.,0.,0.,1.,0.,0.,
0.,0.,0.,0.,0.,0.,
0.,0.,0.,0.,0.,1. },

/*prj_perp[S6P_6]  ,6,*/
{ 1.,0.,0.,0.,0.,0.,
0.,1.,0.,0.,0.,0.,
0.,0.,1.,0.,0.,0.,
0.,0.,0.,1.,0.,0.,
0.,0.,0.,0.,1.,0.,
0.,0.,0.,0.,0.,0. } };




/* The following matrices are the transformation
matrices that may be applied at the associated
boundaries.  *** An addition post permutation may
be needed to preserve a partcular S6 ordering ***
*/

static int S6Rord[24] = { 0,1,2,5,6,7,3,4,14,16,21,23,8,9,10,11,12,13,15,17,18,19,20,22 };

#define S6_Refl_1 0
#define S6_Refl_2 1
#define S6_Refl_3 2
#define S6_Refl_4 3
#define S6_Refl_5 4
#define S6_Refl_6 5
#define S6_Refl_7 6
#define S6_Refl_8 7
#define S6_Refl_9 8
#define S6_Refl_10 9
#define S6_Refl_11 10
#define S6_Refl_12 11
#define S6_Refl_13 12
#define S6_Refl_14 13
#define S6_Refl_15 14
#define S6_Refl_16 15
#define S6_Refl_17 16
#define S6_Refl_18 17
#define S6_Refl_19 18
#define S6_Refl_20 19
#define S6_Refl_21 20
#define S6_Refl_22 21
#define S6_Refl_23 22
#define S6_Refl_24 23

static int S6Refl[24][36] = {

   /*S6Refl_1 B4: { 1 , 2 , 3 , 4 } S6: { 1, 2, 3, 4, 5, 6 ] */
   /*S6_0* identity */
   { 1, 0, 0, 0, 0, 0,
   0, 1, 0, 0, 0, 0,
   0, 0, 1, 0, 0, 0,
   0, 0, 0, 1, 0, 0,
   0, 0, 0, 0, 1, 0,
   0, 0, 0, 0, 0, 1 },

   /*S6Refl_2 B4: { 1 , 2 , 4 , 3 } S6: { 5, 4, 3, 2, 1, 6 } */
   /*S6_1  preserves S6 boundaries 3 and 6*/
{ 0, 0, 0, 0, 1, 0,
0, 0, 0, 1, 0, 0,
0, 0, 1, 0, 0, 0,
0, 1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1 },

/*S6Refl_3 B4: { 1 , 3 , 2 , 4 } S6: { 1, 3, 2, 4, 6, 5 } */
/*S6_2*  preserves boundaries 1 and 4 */
{ 1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 0, 1, 0 },

/*S6Refl_4 B4: { 1 , 3 , 4 , 2 } S6: { 6, 4, 2, 3, 1, 5 } */
/*S6_3  preserves no boundaries */
{ 0, 0, 0, 0, 0, 1,
0, 0, 0, 1, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 1, 0 },

/*S6Refl_5 B4: { 1 , 4 , 2 , 3 } S6: { 5, 3, 4, 2, 6, 1 } */
/*S6_4*/
{ 0, 0, 0, 0, 1, 0,
0, 0, 1, 0, 0, 0,
0, 0, 0, 1, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
1, 0, 0, 0, 0, 0 },

/*S6Refl_6 B4: { 1 , 4 , 3 , 2 } S6: { 6, 2, 4, 3, 5, 1 } */
/*S6_5 preserves boundaries 2 and 5*/
{ 0, 0, 0, 0, 0, 1,
0, 1, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0,
0, 0, 1, 0, 0, 0,
0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0 },

/*S6Refl_7 B4: { 2 , 1 , 3 , 4 } S6: { 2, 1, 3, 5, 4, 6 } */
/*S6_6 preserves boundaries 3 and 6 */
{ 0, 1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 0, 1 },

/*S6Refl_8 B4: { 2 , 1 , 4 , 3 } S6: { 4, 5, 3, 1, 2, 6 } */
/*S6_7 preserves boundaries 3 and 6 */
{ 0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 1, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1 },

/*S6Refl_9 B4: { 2 , 3 , 1 , 4 } S6: { 2, 3, 1, 5, 6, 4 } */
/*S6_8 preserves no boundaries*/
{ 0, 1, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 1, 0, 0 },

/*S6Refl_10 B4: { 2 , 3 , 4 , 1 } S6: { 6, 5, 1, 3, 2, 4 } */
/*S6_9 preserves no boundaries */
{ 0, 0, 0, 0, 0, 1,
0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0 },

/*S6Refl_11 B4: { 2 , 4 , 1 , 3 } S6: { 4, 3, 5, 1, 6, 2 }*/
/*S6_10*/
{ 0, 0, 0, 1, 0, 0,
0, 0, 1, 0, 0, 0,
0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
0, 1, 0, 0, 0, 0 },

/*S6Refl_12 B4: { 2 , 4 , 3 , 1 } S6: { 6, 1, 5, 3, 4, 2 } */
/*S6_11*/
{ 0, 0, 0, 0, 0, 1,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 1, 0, 0, 0,
0, 0, 0, 1, 0, 0,
0, 1, 0, 0, 0, 0 },

/*S6Refl_13 B4: { 3 , 1 , 2 , 4 } S6: { 3, 1, 2, 6, 4, 5  } */
/*S6_12*/
{ 0, 0, 1, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 1, 0 },

/*S6Refl_14 B4: { 3 , 1 , 4 , 2 } S6: { 4, 6, 2, 1, 3, 5 } */
/*S6_13 preserves no boundaries */
{ 0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 0, 1,
0, 1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0,
0, 0, 0, 0, 1, 0 },

/*S6Refl_15 B4: { 3 , 2 , 1 , 4 } S6: { 3, 2, 1, 6, 5, 4 } */
/*S6_14 preserves boundaries 2 and 5*/
{ 0, 0, 1, 0, 0, 0,
0, 1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 0, 1, 0,
0, 0, 0, 1, 0, 0 },

/*S6Refl_16 B4: { 3 , 2 , 4 , 1 } S6: { 5, 6, 1, 2, 3, 4 } */
/*S6_15 */
{ 0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 1,
1, 0, 0, 0, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0,
0, 0, 0, 1, 0, 0 },

/*S6Refl_17 B4: { 3 , 4 , 1 , 2 } S6: { 4, 2, 6, 1, 5, 3 } */
/*S6_16 preserves boundaries 2 and 5 */
{ 0, 0, 0, 1, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 1, 0, 0, 0 },

/*S6Refl_18 B4: { 3 , 4 , 2 , 1 } S6: { 5, 1, 6, 2, 4, 3 } */
/*S6_17*/
{ 0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
0, 1, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0,
0, 0, 1, 0, 0, 0 },

/*S6Refl_19 B4: { 4 , 1 , 2 , 3 } S6: { 3, 5, 4, 6, 2, 1 } */
/*S6_18* preserves no boundaries */
{ 0, 0, 1, 0, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 0, 1,
0, 1, 0, 0, 0, 0,
1, 0, 0, 0, 0, 0 },

/*S6Refl_20 B4: { 4 , 1 , 3 , 2 } S6: ( 2, 6, 4, 5, 3, 1 } */
/*S6_19 preserves no boundaries */
{ 0, 1, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 1, 0, 0, 0,
1, 0, 0, 0, 0, 0 },

/*S6Refl_21 B4: { 4 , 2 , 1 , 3 } S6: { 3, 4, 5, 6, 1, 2 } */
/*S6_20 preserves no boundaries */
{ 0, 0, 1, 0, 0, 0,
0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 1,
1, 0, 0, 0, 0, 0,
0, 1, 0, 0, 0, 0 },

/*S6Refl_22 B4: { 4 , 2 , 3 , 1 } S6: { 1, 6, 5, 4, 3, 2 }*/
/*S6_21 preserves boundaries 1 and 4 */
{ 1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 0, 1, 0,
0, 0, 0, 1, 0, 0,
0, 0, 1, 0, 0, 0,
0, 1, 0, 0, 0, 0 },

/*S6Refl_23 B4: { 4 , 3 , 1 , 2 } S6: {2, 4, 6, 5, 1, 3 } */
/*S6_22  preserves no boundaries */
{ 0, 1, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0 },

/*S6Refl_24 B4: { 4 , 3 , 2 , 1 } S6: { 1, 5, 6, 4, 2, 3 } */
/*S6_23 preserves boundaries 1 and 4*/
{ 1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 1,
0, 0, 0, 1, 0, 0,
0, 1, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0 }

};


#define CS6M_S6Perm_1(v) /* B4:{ 1 , 2 , 3 , 4 }*/  v[0],v[1],v[2],v[3],v[4],v[5]
#define CS6M_S6Perm_2(v) /* B4:{ 1 , 2 , 4 , 3 }*/  v[4],v[3],v[2],v[1],v[0],v[5]
#define CS6M_S6Perm_3(v) /* B4:{ 1 , 3 , 2 , 4 }*/  v[0],v[2],v[1],v[3],v[5],v[4]
#define CS6M_S6Perm_4(v) /* B4:{ 1 , 3 , 4 , 2 }*/  v[5],v[3],v[1],v[2],v[0],v[4]
#define CS6M_S6Perm_5(v) /* B4:{ 1 , 4 , 2 , 3 }*/  v[4],v[2],v[3],v[1],v[5],v[0]
#define CS6M_S6Perm_6(v) /* B4:{ 1 , 4 , 3 , 2 }*/  v[5],v[1],v[3],v[2],v[4],v[0]
#define CS6M_S6Perm_7(v) /* B4:{ 2 , 1 , 3 , 4 }*/  v[1],v[0],v[2],v[4],v[3],v[5]
#define CS6M_S6Perm_8(v) /* B4:{ 2 , 1 , 4 , 3 }*/  v[3],v[4],v[2],v[0],v[1],v[5]
#define CS6M_S6Perm_9(v) /* B4:{ 2 , 3 , 1 , 4 }*/  v[1],v[2],v[0],v[4],v[5],v[3]
#define CS6M_S6Perm_10(v) /* B4:{ 2 , 3 , 4 , 1 }*/  v[5],v[4],v[0],v[2],v[1],v[3]
#define CS6M_S6Perm_11(v) /* B4:{ 2 , 4 , 1 , 3 }*/  v[3],v[2],v[4],v[0],v[5],v[1]
#define CS6M_S6Perm_12(v) /* B4:{ 2 , 4 , 3 , 1 }*/  v[5],v[0],v[4],v[2],v[3],v[1]
#define CS6M_S6Perm_13(v) /* B4:{ 3 , 1 , 2 , 4 }*/  v[2],v[0],v[1],v[5],v[3],v[4]
#define CS6M_S6Perm_14(v) /* B4:{ 3 , 1 , 4 , 2 }*/  v[3],v[5],v[1],v[0],v[2],v[4]
#define CS6M_S6Perm_15(v) /* B4:{ 3 , 2 , 1 , 4 }*/  v[2],v[1],v[0],v[5],v[4],v[3]
#define CS6M_S6Perm_16(v) /* B4:{ 3 , 2 , 4 , 1 }*/  v[4],v[5],v[0],v[1],v[2],v[3]
#define CS6M_S6Perm_17(v) /* B4:{ 3 , 4 , 1 , 2 }*/  v[3],v[1],v[5],v[0],v[4],v[2]
#define CS6M_S6Perm_18(v) /* B4:{ 3 , 4 , 2 , 1 }*/  v[4],v[0],v[5],v[1],v[3],v[2]
#define CS6M_S6Perm_19(v) /* B4:{ 4 , 1 , 2 , 3 }*/  v[2],v[4],v[3],v[5],v[1],v[0]
#define CS6M_S6Perm_20(v) /* B4:{ 4 , 1 , 3 , 2 }*/  v[1],v[5],v[3],v[4],v[2],v[0]
#define CS6M_S6Perm_21(v) /* B4:{ 4 , 2 , 1 , 3 }*/  v[2],v[3],v[4],v[5],v[0],v[1]
#define CS6M_S6Perm_22(v) /* B4:{ 4 , 2 , 3 , 1 }*/  v[0],v[5],v[4],v[3],v[2],v[1]
#define CS6M_S6Perm_23(v) /* B4:{ 4 , 3 , 1 , 2 }*/  v[1],v[3],v[5],v[4],v[0],v[2]
#define CS6M_S6Perm_24(v) /* B4:{ 4 , 3 , 2 , 1 }*/  v[0],v[4],v[5],v[3],v[1],v[2]


static int bdryprsv[6][4] = {
   { 0, 2,21,23 },
{ 0, 5,14,16 },
{ 0, 1,6,7 },
{ 0, 2,21,23 },
{ 0, 5,14,16 },
{ 0, 1,6,7 },
};





static int S6MS[6][36] = {

   /* D6M_1 s3<->s5 {1,2,5,4,3,6} */
   { 1,0,0,0,0,0,
   0,1,0,0,0,0,
   0,0,0,0,1,0,
   0,0,0,1,0,0,
   0,0,1,0,0,0,
   0,0,0,0,0,1 },

   /* D6M_2 s3<->s4 {1,2,4,3,5,6} */

{ 1,0,0,0,0,0,
0,1,0,0,0,0,
0,0,0,1,0,0,
0,0,1,0,0,0,
0,0,0,0,1,0,
0,0,0,0,0,1 },

/* D6M_3 s2<->s4 {1,4,3,2,5,6} */

{ 1,0,0,0,0,0,
0,0,0,1,0,0,
0,0,1,0,0,0,
0,1,0,0,0,0,
0,0,0,0,1,0,
0,0,0,0,0,1 },

/* D6M_4 s2<->s3 {1,3,2,4,5,6} */

{ 1,0,0,0,0,0,
0,0,1,0,0,0,
0,1,0,0,0,0,
0,0,0,1,0,0,
0,0,0,0,1,0,
0,0,0,0,0,1 },

/* D6M_5 s1<->s3 {3,2,1,4,5,6} */

{ 0,0,1,0,0,0,
0,1,0,0,0,0,
1,0,0,0,0,0,
0,0,0,1,0,0,
0,0,0,0,1,0,
0,0,0,0,0,1 },

/* D6M_6  s1<->s2 {2,1,3,4,5,6} */

{ 0,1,0,0,0,0,
1,0,0,0,0,0,
0,0,1,0,0,0,
0,0,0,1,0,0,
0,0,0,0,1,0,
0,0,0,0,0,1 }
};

static int S6MSX[12][36] = {

   /* D6M_1_1 */

   { -1, 0, 0, 0, 0, 0 ,
   1, 1, 0, 0, 0, 0 ,
   1, 0, 0, 0, 1, 0 ,
   -1, 0, 0, 1, 0, 0 ,
   1, 0, 1, 0, 0, 0 ,
   1, 0, 0, 0, 0, 1 },

{ -1, 0, 0, 0, 0, 0 ,
1, 0, 0, 0, 0, 1 ,
1, 0, 1, 0, 0, 0 ,
-1, 0, 0, 1, 0, 0 ,
1, 0, 0, 0, 1, 0 ,
1, 1, 0, 0, 0, 0 },

{ 1, 1, 0, 0, 0, 0 ,
0,-1, 0, 0, 0, 0 ,
0, 1, 0, 1, 0, 0 ,
0, 1, 1, 0, 0, 0 ,
0,-1, 0, 0, 1, 0 ,
0, 1, 0, 0, 0, 1 },

{ 0, 1, 0, 0, 0, 1 ,
0,-1, 0, 0, 0, 0 ,
0, 1, 1, 0, 0, 0 ,
0, 1, 0, 1, 0, 0 ,
0,-1, 0, 0, 1, 0 ,
1, 1, 0, 0, 0, 0 },

{ 1, 0, 1, 0, 0, 0 ,
0, 0, 1, 1, 0, 0 ,
0, 0,-1, 0, 0, 0 ,
0, 1, 1, 0, 0, 0 ,
0, 0, 1, 0, 1, 0 ,
0, 0,-1, 0, 0, 1 },

{ 0, 0, 1, 0, 1, 0 ,
0, 1, 1, 0, 0, 0 ,
0, 0,-1, 0, 0, 0 ,
0, 0, 1, 1, 0, 0 ,
1, 0, 1, 0, 0, 0 ,
0, 0,-1, 0, 0, 1 },

{ 1, 0, 0,-1, 0, 0,
0, 0, 1, 1, 0, 0 ,
0, 1, 0, 1, 0, 0 ,
0, 0, 0,-1, 0, 0 ,
0, 0, 0, 1, 1, 0 ,
0, 0, 0, 1, 0, 1 },

{ 1, 0, 0,-1, 0, 0,
0, 1, 0, 1, 0, 0 ,
0, 0, 1, 1, 0, 0 ,
0, 0, 0,-1, 0, 0 ,
0, 0, 0, 1, 0, 1 ,
0, 0, 0, 1, 1, 0 },

{ 0, 0, 1, 0, 1, 0 ,
0, 1, 0, 0,-1, 0 ,
1, 0, 0, 0, 1, 0 ,
0, 0, 0, 1, 1, 0 ,
0, 0, 0, 0,-1, 0 ,
0, 0, 0, 0, 1, 1 },

{ 1, 0, 0, 0, 1, 0 ,
0, 1, 0, 0,-1, 0 ,
0, 0, 1, 0, 1, 0 ,
0, 0, 0, 0, 1, 1 ,
0, 0, 0, 0,-1, 0 ,
0, 0, 0, 1, 1, 0 },

{ 0, 1, 0, 0, 0, 1 ,
1, 0, 0, 0, 0, 1 ,
0, 0, 1, 0, 0,-1 ,
0, 0, 0, 1, 0, 1 ,
0, 0, 0, 0, 1, 1 ,
0, 0, 0, 0, 0,-1 },

{ 1, 0, 0, 0, 0, 1 ,
0, 1, 0, 0, 0, 1 ,
0, 0, 1, 0, 0,-1 ,
0, 0, 0, 0, 1, 1 ,
0, 0, 0, 1, 0, 1 ,
0, 0, 0, 0, 0,-1 }

};

#define M_1 0
#define M_2 1
#define M_3 2
#define M_4 3
#define M_5 4
#define M_6 5


/* Test if outside S6 region
return 0 if outside */

static int s6test(double s[6]) {
   int retval = 0;
   if (s[0] > 0.) retval |= 0x1;
   if (s[1] > 0.) retval |= 0x2;
   if (s[2] > 0.) retval |= 0x4;
   if (s[3] > 0.) retval |= 0x8;
   if (s[4] > 0.) retval |= 0x10;
   if (s[5] > 0.) retval |= 0x20;

   return retval;
}


/* Compute the dot product of 2 S6 vectors */

static double s6dotprod(double v1[6], double v2[6]) {

   return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3] + v1[4] * v2[4] + v1[5] * v2[5];

}


/* Compute the simple Euclidean distance between 2 S6 vectors */

static double s6eucldistsq(double v1[6], double v2[6]) {

   double vtemp;
   int ii;
   double distsq;

   distsq = 0.;

   for (ii = 0; ii < 6; ii++) {
      vtemp = v1[ii] - v2[ii];
      distsq += vtemp * vtemp;
   }
   return distsq;
}

static double s6eucldist(double v1[6], double v2[6]) {

   return sqrt(s6eucldistsq(v1, v2));

}



#define CS6M_min2(a,b) (((a)<(b))?(a):(b))
#define CS6M_min3(a,b,c) ((CS6M_min2(a,b)<(c))?(CS6M_min2(a,b)):(c))
#define CS6M_min4(a,b,c,d) ((CS6M_min3(a,b,c)<(d))?(CS6M_min3(a,b,c)):(d))
#define CS6M_min5(a,b,c,d,e) ((CS6M_min4(a,b,c,d)<(e))?(CS6M_min4(a,b,c,d)):(e))
#define CS6M_min6(a,b,c,d,e,f) ((CS6M_min5(a,b,c,d,e)<(f))?(CS6M_min5(a,b,c,d,e)):(f))
#define CS6M_min7(a,b,c,d,e,f,g) ((CS6M_min6(a,b,c,d,e,f)<(g))?(CS6M_min6(a,b,c,d,e,f)):(g))
#define CS6M_min8(a,b,c,d,e,f,g,h) ((CS6M_min7(a,b,c,d,e,f,g)<(h))?(CS6M_min7(a,b,c,d,e,f,g)):(h))


/* Macro versions of distances */

#define CS6M_normsq(v) v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]+v[4]*v[4]+v[5]*v[5]
#define CS6M_norm(v) sqrt(CS6M_normsq(v))

#define CS6M_gdistsq(v1,v2) \
( \
(v1[0]-v2[0])*(v1[0]-v2[0])+\
(v1[1]-v2[1])*(v1[1]-v2[1])+\
(v1[2]-v2[2])*(v1[2]-v2[2])+\
(v1[3]-v2[3])*(v1[3]-v2[3])+\
(v1[4]-v2[4])*(v1[4]-v2[4])+\
(v1[5]-v2[5])*(v1[5]-v2[5]))

#define CS6M_gdist(v1,v2) sqrt(CS6M_gdistsq(v1,v2))

#define CS6M_s6eucldistsq(v1,v2) \
( \
(v1[0]-v2[0])*(v1[0]-v2[0])+\
(v1[1]-v2[1])*(v1[1]-v2[1])+\
(v1[2]-v2[2])*(v1[2]-v2[2])+\
(v1[3]-v2[3])*(v1[3]-v2[3])+\
(v1[4]-v2[4])*(v1[4]-v2[4])+\
(v1[5]-v2[5])*(v1[5]-v2[5])\
)

#define CS6M_s6eucldist(v1,v2) sqrt(CS6M_s6eucldistsq(v1,v2))






/*     Compute the best distance between 2 s6 vectors
allowing for permulations of a, b, c and d
*/

#define CS6M_s6eucldistsq_byelem(v11,v12,v13,v14,v15,v16,v21,v22,v23,v24,v25,v26) \
fabs((v11-v21)*(v11-v21)+(v12-v22)*(v12-v22)+(v13-v23)*(v13-v23) + \
(v14-v24)*(v14-v24)+(v15-v25)*(v15-v25)+(v16-v26)*(v16-v26))


#define CS6M_s6prods_byelem(v11,v12,v13,v14,v15,v16,v21,v22,v23,v24,v25,v26) \
2.*(v11*v21+v12*v22+v13*v23+v14*v24+v15*v25+v16*v26)

#define CS6M_s61234distsq(v1,v2) \
CS6M_min3(\
CS6M_min8(\
CS6M_s6eucldistsq_byelem(v1[0],v1[1],v1[2],v1[3],v1[4],v1[5],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[4],v1[3],v1[2],v1[1],v1[0],v1[5],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[0],v1[2],v1[1],v1[3],v1[5],v1[4],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[5],v1[3],v1[1],v1[2],v1[0],v1[4],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[4],v1[2],v1[3],v1[1],v1[5],v1[0],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[5],v1[1],v1[3],v1[2],v1[4],v1[0],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[1],v1[0],v1[2],v1[4],v1[3],v1[5],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[3],v1[4],v1[2],v1[0],v1[1],v1[5],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5])),\
CS6M_min8(\
CS6M_s6eucldistsq_byelem(v1[1],v1[2],v1[0],v1[4],v1[5],v1[3],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[5],v1[4],v1[0],v1[2],v1[1],v1[3],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[3],v1[2],v1[4],v1[0],v1[5],v1[1],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[5],v1[0],v1[4],v1[2],v1[3],v1[1],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[2],v1[0],v1[1],v1[5],v1[3],v1[4],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[3],v1[5],v1[1],v1[0],v1[2],v1[4],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[2],v1[1],v1[0],v1[5],v1[4],v1[3],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[4],v1[5],v1[0],v1[1],v1[2],v1[3],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5])),\
CS6M_min8(\
CS6M_s6eucldistsq_byelem(v1[3],v1[1],v1[5],v1[0],v1[4],v1[2],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[4],v1[0],v1[5],v1[1],v1[3],v1[2],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[2],v1[4],v1[3],v1[5],v1[1],v1[0],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[1],v1[5],v1[3],v1[4],v1[2],v1[0],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[2],v1[3],v1[4],v1[5],v1[0],v1[1],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[0],v1[5],v1[4],v1[3],v1[2],v1[1],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[1],v1[3],v1[5],v1[4],v1[0],v1[2],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5]),\
CS6M_s6eucldistsq_byelem(v1[0],v1[4],v1[5],v1[3],v1[1],v1[2],v2[0],v2[1],v2[2],v2[3],v2[4],v2[5])))

#define CS6M_s61234dist(v1,v2) sqrt(CS6M_s61234distsq(v1,v2))

static double s61234distsq(double v1[6], double v2[6]) {
   double distsq;
   double dwnby;
   double dtrial[24];
   int i;
   distsq = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] + v1[3] * v1[3] + v1[4] * v1[4] + v1[5] * v1[5]
      + v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2] + v2[3] * v2[3] + v2[4] * v2[4] + v2[5] * v2[5];
   dtrial[0] = CS6M_s6prods_byelem(v1[0], v1[1], v1[2], v1[3], v1[4], v1[5], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[1] = CS6M_s6prods_byelem(v1[4], v1[3], v1[2], v1[1], v1[0], v1[5], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[2] = CS6M_s6prods_byelem(v1[0], v1[2], v1[1], v1[3], v1[5], v1[4], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[3] = CS6M_s6prods_byelem(v1[5], v1[3], v1[1], v1[2], v1[0], v1[4], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[4] = CS6M_s6prods_byelem(v1[4], v1[2], v1[3], v1[1], v1[5], v1[0], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[5] = CS6M_s6prods_byelem(v1[5], v1[1], v1[3], v1[2], v1[4], v1[0], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[6] = CS6M_s6prods_byelem(v1[1], v1[0], v1[2], v1[4], v1[3], v1[5], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[7] = CS6M_s6prods_byelem(v1[3], v1[4], v1[2], v1[0], v1[1], v1[5], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[8] = CS6M_s6prods_byelem(v1[1], v1[2], v1[0], v1[4], v1[5], v1[3], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[9] = CS6M_s6prods_byelem(v1[5], v1[4], v1[0], v1[2], v1[1], v1[3], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[10] = CS6M_s6prods_byelem(v1[3], v1[2], v1[4], v1[0], v1[5], v1[1], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[11] = CS6M_s6prods_byelem(v1[5], v1[0], v1[4], v1[2], v1[3], v1[1], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[12] = CS6M_s6prods_byelem(v1[2], v1[0], v1[1], v1[5], v1[3], v1[4], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[13] = CS6M_s6prods_byelem(v1[3], v1[5], v1[1], v1[0], v1[2], v1[4], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[14] = CS6M_s6prods_byelem(v1[2], v1[1], v1[0], v1[5], v1[4], v1[3], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[15] = CS6M_s6prods_byelem(v1[4], v1[5], v1[0], v1[1], v1[2], v1[3], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[16] = CS6M_s6prods_byelem(v1[3], v1[1], v1[5], v1[0], v1[4], v1[2], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[17] = CS6M_s6prods_byelem(v1[4], v1[0], v1[5], v1[1], v1[3], v1[2], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[18] = CS6M_s6prods_byelem(v1[2], v1[4], v1[3], v1[5], v1[1], v1[0], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[19] = CS6M_s6prods_byelem(v1[1], v1[5], v1[3], v1[4], v1[2], v1[0], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[20] = CS6M_s6prods_byelem(v1[2], v1[3], v1[4], v1[5], v1[0], v1[1], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[21] = CS6M_s6prods_byelem(v1[0], v1[5], v1[4], v1[3], v1[2], v1[1], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[22] = CS6M_s6prods_byelem(v1[1], v1[3], v1[5], v1[4], v1[0], v1[2], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dtrial[23] = CS6M_s6prods_byelem(v1[0], v1[4], v1[5], v1[3], v1[1], v1[2], v2[0], v2[1], v2[2], v2[3], v2[4], v2[5]);
   dwnby = dtrial[0];
   for (i = 1; i < 24; i++) if (dtrial[i] > dwnby) dwnby = dtrial[i];
   if (dwnby > distsq) dwnby = distsq;
   return (fabs(distsq - dwnby));
}

#define s61234dist(v1,v2) sqrt(s61234distsq(v1,v2))


static void s6cpyvn(int n, double src[], double dst[]) {
   int i;
   for (i = 0; i < n; i++) {
      dst[i] = src[i];
   }
}


#define CS6M_imv6(v1, m, v2) { \
    int i, j;                  \
    double sum;                \
    for (i = 0; i < 6; i++) {  \
        sum = 0.0;             \
        for(j=0; j < 6; j++) { \
            sum = sum + ((double)m[6*i+j])*v1[j]; \
        }                      \
        v2[i] = sum;           \
    }                          \
}

static void s6imv6(double v1[6], int m[36], double v2[6]) {
   int i, j;
   double sum;
   for (i = 0; i < 6; i++) {
      sum = 0.0;
      for (j = 0; j < 6; j++) {
         sum = sum + ((double)m[6 * i + j])*v1[j];
      }
      v2[i] = sum;
   }
}

#define CS6M_rmv6(v1,m,v2) {   \
    int i, j;                 \
    double sum;               \
    for (i = 0; i < 6; i++) { \
        sum = 0.0;            \
        for(j=0; j < 6; j++) {\
            sum = sum + m[6*i+j]*v1[j]; \
        }                     \
        v2[i] = sum;          \
    }                         \
}

static void s6rmv6(double v1[6], double m[36], double v2[6]) {
   int i, j;
   double sum;
   for (i = 0; i < 6; i++) {
      sum = 0.0;
      for (j = 0; j < 6; j++) {
         sum = sum + m[6 * i + j] * v1[j];
      }
      v2[i] = sum;
   }
}




static void s6twoPminusI(double pg[6], double g[6], double gout[6]) {
   int i;
   for (i = 0; i<6; i++) {
      gout[i] = 2.*pg[i] - g[i];
   }
}


/*     Map a S6 vector onto the boundaries after
applying the 24-way unfolding */


static double s6bddist(double gvec[6], int bdnum) {
   if (bdnum < 0 || bdnum  > 5) return DBL_MAX;
   return gvec[bdnum];
}



/* Get the minimum distance to the boundaries */

static double s6minbddist(double gvec[6]) {
   int ii;
   double dists[6];
   double minbd;
   dists[0] = fabs(gvec[S6P_1]);
   dists[1] = fabs(gvec[S6P_2]);
   dists[2] = fabs(gvec[S6P_3]);
   dists[3] = fabs(gvec[S6P_4]);
   dists[4] = fabs(gvec[S6P_5]);
   dists[5] = fabs(gvec[S6P_6]);

   minbd = dists[0];
   for (ii = 1; ii<6; ii++) {
      minbd = CS6M_min(minbd, dists[ii]);
   }
   return minbd;
}



static void s6bdmaps(double gvec[6],
   double dists[NS6BND],
   int iord[NS6BND],
   double pgs[NS6BND][6],
   double rgs[NS6BND][6],
   double mpgs[NS6BND][6],
   double mvecs[NS6BND][6],
   double maxdist,
   int * ngood) {

   int jj, kk, itemp, igap, idone;

   dists[S6P_1] = fabs(gvec[S6P_1]);
   dists[S6P_2] = fabs(gvec[S6P_2]);
   dists[S6P_3] = fabs(gvec[S6P_3]);
   dists[S6P_4] = fabs(gvec[S6P_4]);
   dists[S6P_5] = fabs(gvec[S6P_5]);
   dists[S6P_6] = fabs(gvec[S6P_6]);


   for (jj = 0; jj < NS6BND; jj++) {
      iord[jj] = jj;
   }


   igap = NS6BND;
   while (igap > 1) {
      igap = igap / 2;
      idone = 0;
      while (!idone) {
         idone = 1;
         for (jj = 0; jj < NS6BND - igap; jj += igap) {
            if (dists[iord[jj]] > dists[iord[jj + igap]]) {
               idone = 0;
               itemp = iord[jj];
               iord[jj] = iord[jj + igap];
               iord[jj + igap] = itemp;
            }
         }
      }
   }

   *ngood = NS6BND;
   for (jj = 0; jj < NS6BND; jj++) {
      for (kk = 0; kk < 6; kk++) {
         pgs[jj][kk] = gvec[kk];
         rgs[jj][kk] = gvec[kk];
      }
      pgs[jj][jj] = 0.;
      rgs[jj][jj] = -rgs[jj][jj];
      CS6M_imv6(pgs[jj], S6MS[jj], mpgs[jj]);
      mpgs[jj][jj] = 0;
      CS6M_imv6(gvec, S6MS[jj], mvecs[jj]);
      if (dists[jj] > maxdist) (*ngood)--;
   }
}

/* revised S6Dist_2bds

Compute lengths of paths from g_lft to g_rgt through bd_up and bd_dwn


*/

static double S6Dist_2bds_rev(double g_lft[6], double g_rgt[6],
   double pg_lft_up[6], double mpg_lft_up[6],
   double pg_rgt_up[6], double mpg_rgt_up[6], int bd_up,
   double pg_rgt_dwn[6], double mpg_rgt_dwn[6],
   double pg_lft_dwn[6], double mpg_lft_dwn[6], int bd_dwn,
   double dist) {

   double bddist_lft_up, bddist_lft_dwn, bddist_right_up, bddist_right_dwn;
   double alpha_lft_up, alpha_rgt_dwn;
   double bdint_lft_up[6], bdint_rgt_dwn[6], mbdint_lft_up[6], mbdint_rgt_dwn[6];
   int ii;
   double gvec_lft_bd_lft_up, gvec_lft_mbd_lft_up, gvec_lft_bd_lft_dwn, gvec_lft_mbd_lft_dwn;
   double gvec_lft_bd_rgt_dwn, gvec_lft_mbd_rgt_dwn;
   double gvec_rgt_bd_lft_up, gvec_rgt_mbd_lft_up, gvec_rgt_bd_lft_dwn, gvec_rgt_mbd_lft_dwn;
   double gvec_rgt_bd_rgt_dwn, gvec_rgt_mbd_rgt_dwn;
   double gvec_lft_bd_lft_upbd_rgt_dwn, gvec_lft_bd_lft_upmbd_rgt_dwn;
   double gvec_lft_mbd_lft_upbd_rgt_dwn, gvec_lft_mbd_lft_upmbd_rgt_dwn;
   double bd_lft_upbd_rgt_dwn, bd_lft_upmbd_rgt_dwn;
   double mbd_lft_upbd_rgt_dwn, mbd_lft_upmbd_rgt_dwn;

   double dist_glft_bd_up, dist_grgt_bd_up;
   double dist_glft_bd_dwn, dist_grgt_bd_dwn;
   double dist_bd_up_bd_dwn, dist_bd_up_mbd_dwn, dist_mbd_up_bd_dwn, dist_mbd_up_mbd_dwn;

   /* If bd_up == bd_dwn:

   We need one mirror point on bd_up

   dividing in proporion to the distance to that boundary

   alpha_lft_up =  g_lft[bd_up]/(g_lft[bd_up]+g_rgt[bd_up]

   if they are different

   We need one mirrow point on bd_up and one on bd_dwn

   alpha_lft_up =  g_lft[bd_up]/(g_lft[bd_up]+g_rgt[bd_up]
   alpha_rgt_dwn =  g_lft[bd_dwn]/(g_lft[bd_dwn]+g_rgt[bd_bdwn]]

   */

   /* distances from g_lft and g_rgt to bd_up and bd_dwn */

   bddist_lft_up = fabs(s6bddist(g_lft, bd_up));
   bddist_right_up = fabs(s6bddist(g_rgt, bd_up));
   if (bd_up != bd_dwn) {
      bddist_lft_dwn = fabs(s6bddist(g_lft, bd_dwn));
      bddist_right_dwn = fabs(s6bddist(g_rgt, bd_dwn));
   }

   /* the portions of a minimal mirror path off each boundary */

   if (bddist_lft_up + bddist_right_up < bddist_lft_up*1.e-10) {
      alpha_lft_up = 0.;
   }
   else {
      alpha_lft_up = bddist_lft_up / (bddist_lft_up + bddist_right_up);
   }
   if (alpha_lft_up > 1.) alpha_lft_up = 1.;
   if (alpha_lft_up < 0.) alpha_lft_up = 0.;
   for (ii = 0; ii < 6; ii++) {
      bdint_lft_up[ii] = pg_lft_up[ii] + alpha_lft_up * (pg_rgt_up[ii] - pg_lft_up[ii]);
   }
   bdint_lft_up[bd_up] = 0.;

   if (bd_up != bd_dwn) {
      if (bddist_lft_dwn + bddist_right_dwn < bddist_lft_dwn*1.e-10) {
         alpha_rgt_dwn = 0.;
      }
      else {
         alpha_rgt_dwn = bddist_lft_dwn / (bddist_lft_dwn + bddist_right_dwn);
      }
      if (alpha_rgt_dwn > 1.) alpha_rgt_dwn = 1.;
      if (alpha_rgt_dwn < 0.) alpha_rgt_dwn = 0.;
      for (ii = 0; ii < 6; ii++) {
         bdint_rgt_dwn[ii] = pg_lft_dwn[ii] + alpha_rgt_dwn * (pg_rgt_dwn[ii] - pg_lft_dwn[ii]);
      }
   }

   CS6M_imv6(bdint_lft_up, S6MS[bd_up], mbdint_lft_up);
   if (bd_up != bd_dwn) {
      CS6M_imv6(bdint_rgt_dwn, S6MS[bd_dwn], mbdint_rgt_dwn);
   }


   gvec_lft_bd_lft_up = s61234dist(g_lft, bdint_lft_up);
   gvec_lft_mbd_lft_up = s61234dist(g_lft, mbdint_lft_up);
   if (gvec_lft_mbd_lft_up < gvec_lft_bd_lft_up)  gvec_lft_bd_lft_up = gvec_lft_mbd_lft_up;

   if (bd_dwn == bd_up) {
      gvec_rgt_bd_lft_up = s61234dist(g_rgt, bdint_lft_up);
      gvec_rgt_mbd_lft_up = s61234dist(g_rgt, mbdint_lft_up);
      if (gvec_rgt_mbd_lft_up <  gvec_rgt_bd_lft_up) gvec_rgt_bd_lft_up = gvec_rgt_mbd_lft_up;
      if (gvec_rgt_mbd_lft_up + gvec_lft_mbd_lft_up < dist) dist = gvec_rgt_mbd_lft_up + gvec_lft_mbd_lft_up;
      return dist;
   }
   else {
      dist_bd_up_bd_dwn = s61234dist(bdint_lft_up, bdint_rgt_dwn);
      dist_bd_up_mbd_dwn = s61234dist(bdint_lft_up, mbdint_rgt_dwn);
      if (dist_bd_up_mbd_dwn < dist_bd_up_bd_dwn) dist_bd_up_bd_dwn = dist_bd_up_mbd_dwn;
      dist_mbd_up_bd_dwn = s61234dist(mbdint_lft_up, bdint_rgt_dwn);
      if (dist_mbd_up_bd_dwn < dist_bd_up_bd_dwn) dist_bd_up_bd_dwn = dist_mbd_up_bd_dwn;
      dist_mbd_up_mbd_dwn = s61234dist(mbdint_lft_up, mbdint_rgt_dwn);
      if (dist_mbd_up_mbd_dwn < dist_bd_up_bd_dwn) dist_bd_up_bd_dwn = dist_mbd_up_mbd_dwn;
      gvec_rgt_bd_rgt_dwn = s61234dist(g_rgt, bdint_rgt_dwn);
      gvec_rgt_mbd_rgt_dwn = s61234dist(g_rgt, mbdint_rgt_dwn);
      if (gvec_rgt_mbd_rgt_dwn < gvec_rgt_bd_rgt_dwn) gvec_rgt_bd_rgt_dwn = gvec_rgt_mbd_rgt_dwn;
      if (gvec_lft_bd_lft_up + dist_bd_up_bd_dwn + gvec_rgt_bd_rgt_dwn < dist) {
         dist = gvec_lft_bd_lft_up + dist_bd_up_bd_dwn + gvec_rgt_bd_rgt_dwn;
      }
      return dist;
   }

   return dist;

}



/*
Compute minimal distance between two Delaunay-reduced
vectors in the Delaunay Cone following the embedding paths
to the 6 boundaries
*/



static double S6Dist_pass(double gvec1[6], double gvec2[6], double dist) {
   double dists1[NS6BND];
   double pgs1[NS6BND][6], rgs1[NS6BND][6], mpgs1[NS6BND][6], mvecs1[NS6BND][6];
   double dists2[NS6BND];
   double pgs2[NS6BND][6], rgs2[NS6BND][6], mpgs2[NS6BND][6], mvecs2[NS6BND][6];
   int iord1[NS6BND], iord2[NS6BND];
   double mindists1;
   double mindists2;
   double disttmp[6], distgather[6], distgather2[6][6];
   int jx1, jx2;
   int ngood1, ngood2;
   double d1[6], d2[6];
   double maxdist;

   maxdist = dist;

   s6bdmaps(gvec1, dists1, iord1, pgs1, rgs1, mpgs1, mvecs1, maxdist, &ngood1);
   s6bdmaps(gvec2, dists2, iord2, pgs2, rgs2, mpgs2, mvecs2, maxdist, &ngood2);

   mindists1 = s6minbddist(gvec1);
   mindists2 = s6minbddist(gvec2);

   if (mindists1 + mindists2 > dist) return dist;

   if (mindists1 + mindists2 < maxdist) {

#pragma omp critical(distminimize)
      {
         for (jx1 = 0; jx1 < 6; jx1++) distgather[jx1] = dist;
      }

#pragma omp parallel for schedule(dynamic)
      for (jx1 = 0; jx1 < ngood1; jx1++) {
         d1[jx1] = dists1[iord1[jx1]];
         disttmp[jx1] = s61234dist(gvec1, pgs1[iord1[jx1]]);
         if (disttmp[jx1] < d1[jx1]) d1[jx1] = disttmp[jx1];
         disttmp[jx1] = s61234dist(gvec1, mpgs1[iord1[jx1]]);
         if (disttmp[jx1] < d1[jx1]) d1[jx1] = disttmp[jx1];
         d2[jx1] = s61234dist(gvec2, pgs1[iord1[jx1]]);
         disttmp[jx1] = s61234dist(gvec2, mpgs1[iord1[jx1]]);
         if (disttmp[jx1] < d2[jx1]) d2[jx1] = disttmp[jx1];
         if (d1[jx1] + d2[jx1] < distgather[jx1]) distgather[jx1] = d1[jx1] + d2[jx1];
      }

#pragma omp flush(dist,distgather)
#pragma omp critical(distminimize)
      {
         for (jx1 = 0; jx1 < 6; jx1++) {
            if (distgather[jx1] < dist) dist = distgather[jx1];
         }
         for (jx1 = 0; jx1 < 6; jx1++) distgather[jx1] = dist;
      }

#pragma omp parallel for schedule(dynamic)
      for (jx2 = 0; jx2 < ngood2; jx2++) {
         d2[jx2] = dists2[iord2[jx2]];
         disttmp[jx2] = s61234dist(gvec2, pgs2[iord2[jx2]]);
         if (disttmp[jx2] < d2[jx2]) d2[jx2] = disttmp[jx2];
         disttmp[jx2] = s61234dist(gvec2, mpgs2[iord2[jx2]]);
         if (disttmp[jx2] < d2[jx2]) d2[jx2] = disttmp[jx2];
         d1[jx2] = s61234dist(gvec1, pgs2[iord2[jx2]]);
         disttmp[jx2] = s61234dist(gvec1, mpgs2[iord2[jx2]]);
         if (disttmp[jx2] < d1[jx2]) d1[jx2] = disttmp[jx2];
         if (d1[jx2] + d2[jx2] < distgather[jx2]) distgather[jx2] = d1[jx2] + d2[jx2];
      }
#pragma omp flush(dist,distgather)
#pragma omp critical(distminimize)
      {
         for (jx1 = 0; jx1 < 6; jx1++) {
            if (distgather[jx1] < dist) dist = distgather[jx1];
         }
         for (jx1 = 0; jx1 < 6; jx1++) distgather[jx1] = dist;
         maxdist = dist;
      }


#pragma omp critical(distminimize)
      {
         for (jx1 = 0; jx1 < 6; jx1++) {
            for (jx2 = 0; jx2 < 6; jx2++) {
               distgather2[jx1][jx2] = dist;
            }
         }
      }

#pragma omp parallel for schedule(dynamic) //collapse(2) removed !!!!
      for (jx1 = 0; jx1 < ngood1; jx1++) {
         for (jx2 = 0; jx2 < ngood2; jx2++) {
            if (dists1[iord1[jx1]] < maxdist && dists2[iord2[jx2]] < maxdist) {
               distgather2[jx1][jx2] = S6Dist_2bds_rev(gvec1, gvec2, pgs1[iord1[jx1]], mpgs1[iord1[jx1]], pgs2[iord1[jx1]], mpgs2[iord1[jx1]], iord1[jx1],
                  pgs2[iord2[jx2]], mpgs2[iord2[jx2]], pgs1[iord2[jx2]], mpgs1[iord2[jx2]], iord2[jx2], distgather2[jx1][jx2]);
            }
         }
      }


#pragma omp flush(dist,distgather2)
#pragma omp critical(distminimize)
      {
         for (jx1 = 0; jx1 < 6; jx1++) {
            for (jx2 = 0; jx2 < 6; jx2++) {
               if (distgather2[jx1][jx2] < dist) dist = distgather2[jx1][jx2];
            }
         }
      }
   }

   return dist;
}



/*
Compute the minimal distance between two S6-reduced
vectors in the S6 Cone following the embedding paths
to the 6 major boundaries
*/



static double CS6Dist(double gvec1[6], double gvec2[6]) {
   double dist, dist1, dist2, distmin;
   dist1 = s6minbddist(gvec1);
   dist2 = s6minbddist(gvec2);
   distmin = CS6M_min(dist1, dist2);
   dist = s61234dist(gvec1, gvec2);
   CSM_report_double("\n  Entered CS6Dist gdist = ", dist, ", ");
   CSM_report_double_vector("gvec1 = ", gvec1, ", ")
      CSM_report_double_vector("gvec2 = ", gvec2, "\n")
      dist = S6Dist_pass(gvec1, gvec2, dist);
   return dist;
}



#endif /*CS6DIST_H */


