/************************************************************************
  atEarthElev.c 	   return value : elevation in radian

  Originally coded                  by C.Otani   1993.2.12
  Modified to use for sisAdSet.c    by T.Dotani
  Modified to use for dp10DyeElv.c  by C.Otani
  Modified                          by C.Otani   1993.12.14
  Modified to atEarthElev           by Y.Ishisaki  94.03.24

	2005/12/04 Y.ISHISAKI	version 2.7
		new algorithm function find_min() to search minimum

	2009/07/24 Y.ISHISAKI	version 3.1
		value CANNOTSEE changed from 120deg -> 200 deg
		stricter xacc condition in while loop of find_min(), adding fabs(y3-y0)

	2009/10/06 Y.ISHISAKI	version 3.2
		bug fix of setting x[1],x[2] in find_min(), pointed out by Bob Wiegand
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "atFunctions.h"
#include "atError.h"

/*#define CANNOTSEE	(120.0*DEG2RAD)*/	/* modified to 200.0 deg in v3.1 */
#define CANNOTSEE	(200.0*DEG2RAD)

/*
 * search minimum of (*f)(x), in the range of xl <= x <= xh
 *
 * NOTE: This function is only valid when (*f)(x) is a convex function.
 *       Searching the minimum at the day/night boundary is the case.
 */
static double				/* ret: x value at minimum */
find_min(
	double (*f)(double x),	/* in: function to minimize */
	double xl,				/* in: lower boundary of x */
	double xh,				/* in: upper boundary of x */
	int maxcall,			/* in: maximum call of (*f)(x) */
	double xacc,			/* in: required accracy for x */
	double *fx				/* out: minimum value of (*f)(x) */
)
{
	double x[4], y[4];
	int imin;
/* flag_hl = 0 if |x3 - x2| > |x1 - x0|, otherwise 1 */
	int flag_hl = 0;

/* start from
   |-----|-----|-----------|  (flag_hl=0)
   x0=xl x1    x2          x3=xh
*/
	x[0] = xl;
	x[2] = (xl + xh) / 2.0;		/* x[2] must be set earlier than x[1] !!! */
	x[1] = (x[0] + x[2]) / 2.0;
	x[3] = xh;

	y[0] = (*f)(x[0]);
	y[1] = (*f)(x[1]);
	y[2] = (*f)(x[2]);
	y[3] = (*f)(x[3]);

	while ( 0 < maxcall-- && xacc < x[3] - x[0] + fabs(y[3] - y[0]) ) {

/* find minimum from x0,x1,x2,x3 */
		imin = 0;
		if ( y[1] < y[imin] ) imin = 1;
		if ( y[2] < y[imin] ) imin = 2;
		if ( y[3] < y[imin] ) imin = 3;

		if ( imin < 2 ) {
			x[3] = x[2];
			y[3] = y[2];
			if ( flag_hl ) {
				x[2] = x[1];
				y[2] = y[1];
				x[1] = (x[0] + x[2]) / 2;
				y[1] = (*f)(x[1]);
/* |-----------|-----|-----|  (flag_hl=1)
   x0          x1    x2    x3
	   to
   |-----|-----|-----|        (flag_hl=0)
   x0    x1    x2    x3 */
			} else {
				x[2] = (x[1] + x[3]) / 2;
				y[2] = (*f)(x[2]);
/* |-----|-----|-----------|  (flag_hl=0)
   x0    x1    x2          x3
	   to
   |-----|--|--|              (flag_hl=1)
   x0    x1 x2 x3 */
			}
		} else {
			x[0] = x[1];
			y[0] = y[1];
			if ( flag_hl ) {
				x[1] = (x[0] + x[2]) / 2;
				y[1] = (*f)(x[1]);
/* |-----------|-----|-----|  (flag_hl=1)
   x0          x1    x2    x3
	   to
               |--|--|-----|  (flag_hl=0)
               x0 x1 x2    x3 */
			} else {
				x[1] = x[2];
				y[1] = y[2];
				x[2] = (x[1] + x[3]) / 2;
				y[2] = (*f)(x[2]);
/* |-----|-----|-----------|  (flag_hl=0)
   x0    x1    x2          x3
	   to
         |-----|-----|-----|  (flag_hl=1)
         x0    x1    x2    x3 */
			}
		}
		flag_hl = ! flag_hl;	/* negate flag_hl, i.e. 0 <-> 1 */
	}

/* find minimum from x0,x1,x2,x3 */
	imin = 0;
	if ( y[1] < y[imin] ) imin = 1;
	if ( y[2] < y[imin] ) imin = 2;
	if ( y[3] < y[imin] ) imin = 3;

	*fx = y[imin];
/*	printf("\
find_min: x=%.6f, f(x)=%.6f, acc=%.4e maxcall=%d\n",
		x[imin], *fx, x[3]-x[0], maxcall);*/

	return x[imin];
}

#if 0
static double
func(double x)
{
	AtVect vect, nvect;
	double fact;
	vect[0] = eradw * cos(x) - satw[0];
	vect[1] = eradw * sin(x);
	vect[2] = -satw[2];
	atNormVect(vect, nvect);
	return - atScalProd(sisw, nvect);
}
#endif

static double funcA;	/* vSat[0]* nvTgt[0] + vSat[2]*nvTgt[2] */
static double funcB;	/* erad * nvTgt[0] */
static double funcC;	/* erad * nvTgt[1] */
static double funcD;	/* erad**2 + vSat[0]**2 + vSat[2]**2 */
static double funcE;	/* 2 * erad * vSat[0] */

static double
fastfunc(double x)
{
	double cosx = cos(x);
	double sinx = sin(x);
	return (funcA - funcB*cosx - funcC*sinx) / sqrt(funcD - funcE*cosx);
}

/*
 * numcal2() search minimum of fastfunc(),
 * which is invert of the scaler product of
 * normalized vector of the field of view (nvTgt), and
 * normalized vector to the border of day/night on the Earth surface
 * from the satellite. The coordinate of center of the Earth is
 * E (-vSat[0], 0.0==vSat[1], -vSat[2]), and the vector
 * EP = (erad*cos(x), erad*sin(x), 0.0) = erad (cos(x), sin(x), 0.0)
 * OP = (erad*cos(x) - vSat[0], erad*sin(x), 0.0),
 * where - bound <= x <= + bound.
 * At the boundary x = b, (OP . EP) = 0.0, therefore,
 *	erad cos(b) cos(b) - vSat[0] cos(b) + erad sin(b) sin(b) = 0
 * 	cos(b) = erad / vSat[0]
 */
static double
numcal2(AtVect vSat, AtVect nvTgt, double erad)
{
	static int maxcall = 100;
	static double xacc = 1e-9;
	double bound, x, y, ang;

	funcA = vSat[0]*nvTgt[0] + vSat[2]*nvTgt[2];
	funcB = erad * nvTgt[0];
	funcC = erad * nvTgt[1];
	funcD = erad*erad + vSat[0]*vSat[0] + vSat[2]*vSat[2];
	funcE = 2 * erad * vSat[0];
	if ( vSat[0] < erad ) return CANNOTSEE;
	bound = acos(erad / vSat[0]);
/*	printf("bound=%.6f\n", bound*RAD2DEG);*/

	x = find_min(fastfunc, -bound, +bound, maxcall, xacc, &y);

	if ( y < -1.0 ) {
		y = -1.0;
	} else if ( 1.0 < y ) {
		y = +1.0;
	}
	ang = acos(-y);

	return ang;
}

int
atEarthElev(
	AtVect vSat,	/*  (in) sidereal   */
	AtVect nvTgt,	/*  (in) sidereal   */
	AtVect vSun,	/*  (in) sidereal   */
	int *earth_occult,	/* (out) earth occultation */
	double elevation[3]	/* (out) from earth, DE, NE */
)
{
	static double erad = EARTH_RADIUS;

	AtVect vEarth, vSat2, vTgt2;
	AtRotMat RM;
	double elv, angEarth, dAng, xDist, cross, dot, satDist;
	double xyradius;

	/*   elevation angle   */
	ATInvVect(vSat, vEarth);
	satDist = ATNorm(vEarth);
	angEarth = asin(erad/satDist);
	atAngDistance(nvTgt, vEarth, &xDist);
	elevation[0] = elv = xDist - angEarth;

	/*   coordinate transformation   */
	if ( atSetRotMatZX(vSun, vSat, RM) ) {
		static AtVect xaxis = { 1.0, 0.0, 0.0 };
		static AtVect yaxis = { 0.0, 1.0, 0.0 };
		static AtVect zaxis = { 0.0, 0.0, 1.0 };
		if ( 0.0 == satDist || 0.0 == ATNorm(vSun) ) {
			return NULL_VECTOR;
		}
		/* vSun // vSat */
		if ( atSetRotMatZX(vSun, xaxis, RM) ) {
			if ( atSetRotMatZX(vSun, yaxis, RM) ) {
				if ( atSetRotMatZX(vSun, zaxis, RM) ) {
					return NULL_VECTOR;
				}
			}
		}
	}
	ATRotVect(RM, vSat, vSat2);
	ATRotVect(RM, nvTgt, vTgt2);

/* after the rotation, the Sun is at the +Z-axis direction
   from the satellite, the Earth is on the X-Z plane where x <= 0.0.
   i.e. 0.0 <= vSat2[0] && 0.0 == vSat2[1] because vSat2 is the vector
   to the satellite from center of the Earth.
   The Earth is lit by the Sun for +Z direction, night at -Z direction. */

	/* If xyradius <= eradius, the result is very simple. */
	xyradius = vSat2[0];       /*  always vSat2[1]==0.0 */
	if ( xyradius <= erad ) {
		if ( vSat2[2] <= 0.0  ) {
			/*  night earth */
			elevation[1] = CANNOTSEE;
			elevation[2] = elv;
		} else {
			/* bright earth */
			elevation[1] = elv;
			elevation[2] = CANNOTSEE;
		}
	} else if ( elv <= 0.0 ) {
		/* satellite is looking at the Earth */
/*
   ->
   p : vector to the cross point on the Earth surface [ unknown ]

   ->
   g : normalized vector to the target from the satellite [ vTgt2 ]

   ->
   e : vector to the center of the Earth from the satellite [ - vSat2 ]

     ->          ->           ->  ->    ->
   | g | = 1,  | p | = erad,  e + p = k g   (1)

              ->    -> 2  ->  ->     ->  ->
   (1) times  e : | e | + e . p  = k e . g  (2)

              ->    -> 2  ->  ->     ->  ->
   (1) times  p : | p | + e . p  = k p . g  (3)

                    -> 2    -> 2                                  -> 2
   (2) + (3)    : | e | + | p | + 2 e . p = k g . (e + p) = k^2 | g |

                            ->  ->    -> 2    -> 2
   using (2)    : k^2 - 2 k e . g + | e | - | p |  = 0

                      ->  ->               ->  -> 2     -> 2    -> 2
   solving k    : k = e . g   +/-  sqrt[ ( e . g )  + | p | - | e |  ]

   choosing smaller :  k = dot - sqrt(dot*dot + erad*erad - satDist*satDist)

                  ->    ->  ->
   therefore    : p = k g - e = k * vTgt2 + vSat2

*/
		dot = ATScalProd(vEarth, nvTgt);	/* dot is invariant for rotation */
		cross = vSat2[2];
		cross += vTgt2[2] * (dot - sqrt(erad*erad-satDist*satDist+dot*dot));
		if ( 0.0 <= cross ) {
			/* looking at the bright earth */
			vSat2[2] = -vSat2[2];
			vTgt2[2] = -vTgt2[2];
			dAng = numcal2(vSat2, vTgt2, erad);
			elevation[1] = -dAng;
			elevation[2] = dAng;
		} else {
			/* looking at the dark earth */
			dAng = numcal2(vSat2, vTgt2, erad);
			elevation[1] = dAng;
			elevation[2] = -dAng;
		}
	} else {
		/* looking at the space */
/*
   Point D represents the vertical position from center of the Earth
   in the line of sight to the target [vTgt2], and the vector vD
   represents the vector to the Point D from center of the Earth.
   The distance to D from the satellite is
           ->  ->
   dot = | e . g | = satDist | cos(theta) |

   The point D may be an opossite side of the target vector,
   The following formulae can evaluate day/night side correctly.
*/
/*		dot = fabs(ATScalProd(vSat2,vTgt2));
		vD[0] = vSat2[0] + dot*vTgt2[0];
		vD[1] =            dot*vTgt2[1];
		vD[2] = vSat2[2] + dot*vTgt2[2];
		cross = ( vD[2] / ATNorm(vD) ) * erad;*/

		dot = fabs(ATScalProd(vSat, nvTgt));/* dot is invariant for rotation */
		cross = vSat2[2] + dot * vTgt2[2];	/* scaling is not required */
		if ( cross >= 0.0 ) {
			/* near the bright earth */
			vSat2[2] = -vSat2[2];
			vTgt2[2] = -vTgt2[2];
			dAng = numcal2(vSat2, vTgt2, erad);
			elevation[1] = elv;
			elevation[2] = dAng;
		} else {
			/* near the dark earth */
			dAng = numcal2(vSat2, vTgt2, erad);
			elevation[1] = dAng;
			elevation[2] = elv;
		}
	}
	if ( 0 < elv ) {
		*earth_occult = 0;	/* sky */
	} else if ( 0 < elevation[1] ) {
		*earth_occult = 1;	/* night earth */
	} else {
		*earth_occult = 2;	/* bright earth */
	}
	return 0;
}

#if 0	/* comment in Japanese */

Date: Fri, 11 Nov 2005 02:51:57 +0900 (JST)
From: Yoshitaka ISHISAKI <ishisaki@phys.metro-u.ac.jp>
Subject: [suzakuhelp] Re: rev0.3 criteria
To: suzakuhelp@astro.isas.jaxa.jp

堂谷さま、みなさま

>> On Fri, 11 Nov 2005 00:02:25 +0900,
>> 堂谷 忠靖 <dotani@astro.isas.jaxa.jp> said:
>> もう一度実際のehkファイルを見てみたのですが、両者は
>> 必ずしも一致しないようです。極端な例ですが、添付した
>> ファイルをご覧下さい。赤がELV、緑がNTE_ELV、
>> 青がDYE_ELVです。
>>
> どういうアルゴリズムでNTE_ELV等が計算されているか、
> すっかり忘れてしまいましたが、ひょっとしたらELV<0のときは、
> 地球の縁を無視して、つねに地球上の昼夜の境目からの
> 離角を計算するようになっているかもしれません。

計算は atFunctions の atEarthElev.c で行なっていますが、
確かにそのようになっています。従って、夜地球/昼地球の
切り出しには NTE_ELV/DYE_ELV だけでなく、ELV も併用すべきです。

判断の目安とするため、それぞれの場合に対する設定値の定義を
書いておきます。

1. 衛星から見える地球が完全に昼の時
	DYE_ELEV = ELV
	NTE_ELEV = 120 deg

2. 衛星から見える地球が完全に夜の時
	DYE_ELEV = 120 deg
	NTE_ELEV = ELV

3. 衛星から見える地球が昼/夜に分割されている時

   3-a  衛星が昼地球を向いている時
	DYE_ELEV = - |視野方向と地球の昼/夜の境目からの角度|
	NTE_ELEV = + |視野方向と地球の昼/夜の境目からの角度|

   3-b  衛星が夜地球を向いている時
	DYE_ELEV = + |視野方向と地球の昼/夜の境目からの角度|
	NTE_ELEV = - |視野方向と地球の昼/夜の境目からの角度|

   3-c  衛星が空を向いていて、一番近い地球のへりが昼の時
	DYE_ELEV = ELEV
	NTE_ELEV = |視野方向と地球の昼/夜の境目からの角度|

   3-d  衛星が空を向いていて、一番近い地球のへりが夜の時
	DYE_ELEV = |視野方向と地球の昼/夜の境目からの角度|
	NTE_ELEV = ELEV

---
石崎 欣尚 (ISHISAKI Yoshitaka)
首都大学東京/東京都立大学 理学研究科  TEL: 0426-77-2485  FAX: 0426-77-2483
E-Mail: ishisaki@phys.metro-u.ac.jp

#endif
