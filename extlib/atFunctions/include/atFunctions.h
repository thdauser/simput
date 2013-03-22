/*************************************************************************
	type definitions of "At" library

  1993/01/26 N.Kawai		v1.5
  1996/09/20 Y.ISHISAKI		v1.8
  2004/03/10 Y.ISHISAKI		v2.2
	add AtTimeD structure & related routines, mission time
  2004/07/06 Y.ISHISAKI		v2.3
	add declaration of AtElement structure
  2005/08/06 Y.ISHISAKI		v2.5
	EARTH_RADIUS, EARTH_E2 are modified to use GRS80
	add declarations of atGeodeticToGeographic(), atGeographicToGeodetic()
  2005/08/06 Y.ISHISAKI		v2.5
	add declarations of atInterpolateQuat(), atInterpolateEuler()
  2005/08/14 Y.ISHISAKI		v2.5
  	add declaration of atMissionTimeResetTable()
  2005/08/17 Y.ISHISAKI		v2.5
	add declaration of atQuatToEuler(), atEulerToQuat()
  2005/10/09 Y.ISHISAKI		v2.5
	add declarations of atGeodesic(), atHXDBrazil(), atPathNumUSC()
  2005/12/04 Y.ISHISAKI		v2.7
  	add atSubVect(), atMulVect(), atDivVect()
	add fast macros, ATCopyVect, ATInvVect,
		ATAddVect, ATSubVect(), ATMulVect, ATDivVect,
		ATMulAddVect, ATVectProd, ATScalProd, ATNorm,
		ATRotVect, ATInvRotMat, ATRMProd
	add atRigidityD(), now atRigSet(), atRigidity() is obsolete
  2006/02/08 Y.ISHISAKI		v2.8
	add declaration of atSAA()
	add comments that longitude & latitude are in geodetic coordinate,
	for atBrazil(), atSISBrazil(), atSTTBrazil(), atHXDBrazil()
	fix comment for atRigidityD(), atRigidity(), geographic -> geodetic
	define SAA_NOMINAL,SAA_ASCA_SIS,SAA_ASCA_STT,SAA_SUZAKU_HXD
  2009/07/24 Y.ISHISAKI		v3.1
	add new typedef of AtElement3
	add new functions, atSetElement3(), atPathNum3(), atPathNumUSC3(),
		atSatPos3(), atElementTime3()
**************************************************************************/

#ifndef _ATFUNCTIONS_H_
#define _ATFUNCTIONS_H_

/* fast macros, since atFunctions-2.7 */

/* same as atCopyVect() */
#define ATCopyVect(X,Y)	\
((Y)[0]=(X)[0],(Y)[1]=(X)[1],(Y)[2]=(X)[2],0)
/* same as atInvVect() */
#define ATInvVect(X,Y) \
((Y)[0]=-(X)[0],(Y)[1]=-(X)[1],(Y)[2]=-(X)[2],0)
/* same as atAddVect() */
#define ATAddVect(X,Y,Z) \
((Z)[0]=(X)[0]+(Y)[0],(Z)[1]=(X)[1]+(Y)[1],(Z)[2]=(X)[2]+(Y)[2],0)
/* same as atSubVect() */
#define ATSubVect(X,Y,Z) \
((Z)[0]=(X)[0]-(Y)[0],(Z)[1]=(X)[1]-(Y)[1],(Z)[2]=(X)[2]-(Y)[2],0)
/* same as atMulVect() */
#define ATMulVect(F,X,Z) \
((Z)[0]=(F)*(X)[0],(Z)[1]=(F)*(X)[1],(Z)[2]=(F)*(X)[2],0)
/* same as atDivVect() */
#define ATDivVect(D,X,Z) \
(0.0==(D)?-1:((Z)[0]=(X)[0]/(D),(Z)[1]=(X)[1]/(D),(Z)[2]=(X)[2]/(D),0))
/* same as atMulAddVect() */
#define ATMulAddVect(F,X,G,Y,Z)	\
((Z)[0]=(F)*(X)[0]+(G)*(Y)[0],\
 (Z)[1]=(F)*(X)[1]+(G)*(Y)[1],\
 (Z)[2]=(F)*(X)[2]+(G)*(Y)[2],\
0)
/* same as atVectProd() */
#define ATVectProd(X,Y,Z) \
((Z)[0]=(X)[1]*(Y)[2]-(X)[2]*(Y)[1],\
 (Z)[1]=(X)[2]*(Y)[0]-(X)[0]*(Y)[2],\
 (Z)[2]=(X)[0]*(Y)[1]-(X)[1]*(Y)[0],\
0)
/* same as atScalProd() */
#define ATScalProd(X,Y) \
((X)[0]*(Y)[0]+(X)[1]*(Y)[1]+(X)[2]*(Y)[2])
/* same as atNorm() */
#define ATNorm(V) \
(sqrt(ATScalProd(V,V)))
/* same as atRotVect() */
#define ATRotVect(M,X,Y) \
((Y)[0]=(M)[0][0]*(X)[0]+(M)[0][1]*(X)[1]+(M)[0][2]*(X)[2],\
 (Y)[1]=(M)[1][0]*(X)[0]+(M)[1][1]*(X)[1]+(M)[1][2]*(X)[2],\
 (Y)[2]=(M)[2][0]*(X)[0]+(M)[2][1]*(X)[1]+(M)[2][2]*(X)[2],\
0)
/* same as atInvRotMat() */
#define ATInvRotMat(RM1,RM2) \
((RM2)[0][0]=(RM1)[0][0],(RM2)[0][1]=(RM1)[1][0],(RM2)[0][2]=(RM1)[2][0],\
 (RM2)[1][0]=(RM1)[0][1],(RM2)[1][1]=(RM1)[1][1],(RM2)[1][2]=(RM1)[2][1],\
 (RM2)[2][0]=(RM1)[0][2],(RM2)[2][1]=(RM1)[1][2],(RM2)[2][2]=(RM1)[2][2],\
0)
/* same as atRMProd() */
#define ATRMProd(M0,M1,M2) \
((M2)[0][0]=(M1)[0][0]*(M0)[0][0]+(M1)[0][1]*(M0)[1][0]+(M1)[0][2]*(M0)[2][0],\
 (M2)[0][1]=(M1)[0][0]*(M0)[0][1]+(M1)[0][1]*(M0)[1][1]+(M1)[0][2]*(M0)[2][1],\
 (M2)[0][2]=(M1)[0][0]*(M0)[0][2]+(M1)[0][1]*(M0)[1][2]+(M1)[0][2]*(M0)[2][2],\
 (M2)[1][0]=(M1)[1][0]*(M0)[0][0]+(M1)[1][1]*(M0)[1][0]+(M1)[1][2]*(M0)[2][0],\
 (M2)[1][1]=(M1)[1][0]*(M0)[0][1]+(M1)[1][1]*(M0)[1][1]+(M1)[1][2]*(M0)[2][1],\
 (M2)[1][2]=(M1)[1][0]*(M0)[0][2]+(M1)[1][1]*(M0)[1][2]+(M1)[1][2]*(M0)[2][2],\
 (M2)[2][0]=(M1)[2][0]*(M0)[0][0]+(M1)[2][1]*(M0)[1][0]+(M1)[2][2]*(M0)[2][0],\
 (M2)[2][1]=(M1)[2][0]*(M0)[0][1]+(M1)[2][1]*(M0)[1][1]+(M1)[2][2]*(M0)[2][1],\
 (M2)[2][2]=(M1)[2][0]*(M0)[0][2]+(M1)[2][1]*(M0)[1][2]+(M1)[2][2]*(M0)[2][2],\
0)

/*
 * 	Three-dimensional Cartesian Vector.
 */
typedef double AtVect[3];

/*
 * 	3 x 3 Rotation Matrix.
 */
typedef double AtRotMat[3][3];

/*
 * 	Euler angle notation of rotation. (usually z-y-z)
 */
typedef struct {
	double phi; 	/* First Euler Angle (radian)*/
	double theta; 	/* Second Euler Angle (radian)*/
	double psi;	/* Third Euler Angle (radian)*/
} AtEulerAng;

/*
 * 	Three-dimensional Vector in Polar coordinate.
 */
typedef struct {
	double r;	/* radial component */
	double lon; 	/* longitude or R.A. (radian)*/
	double lat;	/* latitude or declination (radian)*/
} AtPolarVect;

typedef struct {
	unsigned int hour;	/* hour part  */
	unsigned int min; 	/* minuit */
	double sec;		/* second */
} AtRightAscension;

typedef struct {
	int sign;              /* +1: north or -1:south */
	unsigned int deg;      /* degree part  */
	unsigned int min;      /* minuit */
	double sec;           /* second */
} AtDeclination;

typedef struct {
	AtRightAscension ra;	/* right ascension in hhmmss  */
	AtDeclination dec; 	/* declination in ddmmss */
} AtPolarVect60;

/*
 * 	Quaternion notation of rotation (Q-parameters).
 */
typedef double AtQuat[4];

/*
 *	Time in year, month, day, hour, minute, sec
 */
typedef struct {
	int yr;		/* year  */
	int mo; 	/* month */
	int dy; 	/* day */
	int hr; 	/* hour */
	int mn; 	/* minute */
	int sc; 	/* sec */
	float ms; 	/* ms */
} AtTime;

typedef struct {
	int yr;		/* year  */
	int mo; 	/* month */
	int dy; 	/* day */
	int hr; 	/* hour */
	int mn; 	/* minute */
	int sc; 	/* sec */
	double ss; 	/* sub-sec in unit of sec */
} AtTimeD;

/*
 *	Orbital Elements (angles are in radian), OBSOLETE & INTERNAL USE ONLY
 */
typedef struct {
    AtTime itz;	/* Epoch of Elements in year, month, day, hour, min, sec */
    double semiax;	/* Semi-major Axis in km */
    double eccent;	/* Eccentricity */
    double aincln;	/* Inclination */
    double ragome;	/* Right Ascension of Ascending Node */
    double smaome;	/* Argument of Perigee (angle from ascending node) */
    double omean0;	/* Mean Anomaly at the Epoch */
    double adot;	/* Time Derivative of "semiax" in km/day */
    double eccdot;      /* Time Derivative of "eccent" in /day (v1.6)*/
    double aindot;      /* Time Derivative of "aincln" in radian/day (v1.6)*/
    double ragdot;	/* Time Derivative of "ragome" in radian/day */
    double smodot;	/* Time Derivative of "smaome" in radian/day */
    double znbar;	/* Time Derivative of Mean Anomaly in radian/min,
    			i.e. Mean Motion */
    double znbadt;	/* 2nd Derivative of Mean Anomaly in radian/min/day */
    double mjdz;	/* Modified Julian Day of the Epoch */
    double perige;	/* Perigee */
    double apoge;	/* Apogee */
} AtElement;

/*
 *	Orbital Elements (angles are in radian)
 */
typedef struct {
    AtTimeD itz;	/* Epoch of Elements in yr, mo, dy, hr, mn, sc, ss */
    AtTimeD itn;	/* next Epoch of Elements */
    double mjdz;	/* Modified Julian Day of the Epoch */
    double mjdn;	/* next Epoch in MJD */
    double semiax;	/* Semi-major Axis in km */
    double eccent;	/* Eccentricity */
    double aincln;	/* Inclination */
    double ragome;	/* Right Ascension of Ascending Node */
    double smaome;	/* Argument of Perigee (angle from ascending node) */
    double omean0;	/* Mean Anomaly at the Epoch */
    double adot;	/* Time Derivative of "semiax" in km/day */
    double eccdot;      /* Time Derivative of "eccent" in /day (v1.6)*/
    double aindot;      /* Time Derivative of "aincln" in radian/day (v1.6)*/
    double ragdot;	/* Time Derivative of "ragome" in radian/day */
    double smodot;	/* Time Derivative of "smaome" in radian/day */
    double znbar;	/* Time Derivative of Mean Anomaly in radian/min,
    			i.e. Mean Motion */
    double znbadt;	/* 2nd Derivative of Mean Anomaly in radian/min/day */
    double perige;	/* Perigee */
    double apoge;	/* Apogee */
} AtElement3;

/*
 *	Cut-off rigidity data
 */
typedef struct {
	char *filename;
	int nx;
	int ny;
	double altitude;	/* altitude [km] from center of the Earth */
	double *slon;		/* slon[nx] in geocentric coordinates */
	double *slat;		/* slat[ny] in geocentric coordinates */
	double *coeff;		/* coeff[nx][ny] */
	double *data;		/* data[nx][ny] */
	double work[1];		/* work[5*nx] */
} AtRigData2;

/**************************************************************************/
/*									*/
/*	Global Variable (COMMON variables in FORTRAN)			*/
/*									*/
/**************************************************************************/

extern AtElement atElement;

/**************************************************************************/
/*									*/
/*	Basic constants							*/
/*									*/
/**************************************************************************/

#define PI		3.1415926535897932385
#define TWO_PI		6.283185307179586477
#define DEG2RAD		0.017453292519943295769
#define RAD2DEG		57.295779513082320877

/* GRS80, http://dgfi2.dgfi.badw-muenchen.de/geodis/REFS/grs80.html */
#define EARTH_RADIUS	6378.137
#define EARTH_E2	0.00669438002290

#define MOON_RADIUS	1738.
#define AU		149597870.
#define EPS		1.e-12
#define MJD_B1950	33281.923
#define MJD_J2000	51544.500

#define SAA_NOMINAL	0
#define SAA_ASCA_SIS	1
#define SAA_ASCA_STT	2
#define SAA_SUZAKU_HXD	3

/**************************************************************************/
/*									*/
/*	Basic unit conversion						*/
/*									*/
/**************************************************************************/

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * covert Right Ascension in hr,min,s to radian.
 */
double atRAToRadian(		/* right ascension in radian  */
	AtRightAscension ra);	/* right ascension in hhmmss  */

/*
 * covert Declination in deg,min,s to radian.
 */
double atDecToRadian(		/* declination in radian  */
	AtDeclination dec); 	/* declination in ddmmss */

/*
 * covert R.A. and Dec in hr/deg,min,s to Vector.
 */
int atPol60ToVect(
	    AtPolarVect60 *p, 	/* input: polar vector in hh/deg mm ss */
	    AtVect x);		/* output: vector */

/*
 * covert Vector to R.A. and Dec in hr/deg,min,s.
 */
int atVectToPol60(
	    AtVect x,		/* input: vector */
	    AtPolarVect60 *p); 	/* output: polar vector in hh/deg mm ss */

/*
 * covert Vector to R.A. and Dec in degree.
 */
int atVectToPolDeg(
	    AtVect x,		/* input: vector */
	    double *r, 		/* vector length */
	    double *alpha, 	/* R.A. in degree */
	    double *delta); 	/* Dec. in degree */

/*
 * covert R.A. and Dec in degree to Vector
 */
int atPolDegToVect(
	    double r, 		/* vector length */
	    double alpha, 	/* R.A. in degree */
	    double delta, 	/* Dec. in degree */
	    AtVect x);		/* output: vector */

/************************************************************************/
/*									*/
/*	Basic coordinate conversion subroutines				*/
/*									*/
/************************************************************************/

/*
 * ROTATE A VECTOR WITH ROTATION MATRIX.
 */
int atRotVect(
	AtRotMat rm,	/* input: rotation matrix */
	AtVect x,	/* input: vector */
	AtVect y);	/* output: vector */

/*
 * CALC ROTATION MATRIX FOR AXIS AND ROLL ANGLE.
 */
int atSetRotMat(
	AtVect axis,	/* input: axis of rotation, should be non zero */
	double roll,	/* input: roll angle around axis (radian) */
	AtRotMat rm);	/* output: rotation matrix */

/*
 * ROTATION MATRIX defined by New Z-axis and a vector in (+)X-Z half plane.
 */
int atSetRotMatZX(
	AtVect zAxis,	/* input: vector defining new z-axis */
	AtVect xAxis,	/* input: vector in new +X-Z half plane */
	AtRotMat rm);	/* output: rotation matrix */

/*
 * CALC INVERSE ROTATION MATRIX.
 */
int atInvRotMat(
	AtRotMat rm,	/* input: rotation matrix */
	AtRotMat rm2);	/* output: inversed rotation matrix */

/*
 * NORMALIZE A VECTOR.
 */
int atNormVect(
	AtVect x,	/* input: vector */
	AtVect y);	/* output: normalized vector*/

/*
 * Norm (Absolute size) of A VECTOR.
 */
double atNorm(	/* output: normalized vector*/
	AtVect x);	/* input: vector */

/*
 * product of two rotation matrices rm2 = rm1 rm0
 */

int atRMProd(
	AtRotMat rm0,	/* input: rotation matrix to be multiplied*/
	AtRotMat rm1,	/* input: rotation matrix to multiply*/
	AtRotMat rm2);	/* output: product rm2 = rm1 rm0 */

/*
 * Checking consistency (unitarity) of a rotation matrix
 */

int atRMCheck(
	AtRotMat rm);	/* input/output: rotation matrix to be checked*/

/*
 * FIND EULER ANGLES FOR A ROTATION MATRIX.
 */
int atRMToEuler(
	AtRotMat rm,		/* input: rotation matrix */
	AtEulerAng *ea);	/* output: z-y-z Euler Angle (radian) */

/*
 * FIND ROTATION MATRIX FOR A EULER ANGLE SET.
 */
int atEulerToRM(
	AtEulerAng *ea,		/* input: z-y-z Euler Angle (radian) */
	AtRotMat rm);		/* output: rotation matrix */

/*
 * convert a Rotation Matrix to Quaternion.
 */
int atRMToQuat(
	AtRotMat rm,	/* input: rotation matrix */
	AtQuat q);	/* output: quaternion */

/*
 * FIND ROTATION MATRIX FOR A EULER ANGLE SET.
 */
int atQuatToRM(
	AtQuat q,		/* input: quaternion */
	AtRotMat rm);		/* output: rotation matrix */

/*
 * product of two quaternion: q2 = q0 q1
 * 	(in matrix representation: rm(q2) = rm(q1) rm(q0).
 *	note the inverse order of quaternion as compared to matrix )
 */

int atQuatProd(
	AtQuat q0,	/* input: quaternion to be multiplied*/
	AtQuat q1,	/* input: quaternion to multiply*/
	AtQuat q2);	/* output: product */

/************************************************************************/
/*									*/
/*	Attitude and Orbit functions Library based on CSUB.FORT77	*/
/*									*/
/************************************************************************/

/* HENKAN */

/*
 * copy vectors (z = x)
 */
int atCopyVect(
	AtVect x,	/* in: vector */
	AtVect z);	/* out: copy vector of x*/

/*
 * calc inverse vector  (RVSVCT)
 */
int atInvVect(
	AtVect x,	/* in: vector */
	AtVect y);	/* out: reversed vector */

/*
 * add vectors (z = x + y)
 */
int atAddVect(
	AtVect x,	/* in: vector */
	AtVect y,	/* in: vector */
	AtVect z);	/* out: sum vector of x and y */

/*
 * subtract vectors (z = x - y)
 */
int atSubVect(
	AtVect x,	/* in: vector */
	AtVect y,	/* in: vector */
	AtVect z);	/* out: vector of (x - y) */

/*
 * multiply scaler for vector (z = f*x)
 */
int atMulVect(
	double f,	/* in: multiplicand for x */
	AtVect x,	/* in: vector */
	AtVect z);	/* out: answer*/

/*
 * divide scaler for vector (z = (1/d)*x)
 */
int atDivVect(
	double d,	/* in: divisor for x */
	AtVect x,	/* in: vector */
	AtVect z);	/* out: answer*/

/*
 * multiply and add vectors (z = f*x + g*y)
 */
int atMulAddVect(
	double f,	/* in: multiplicand for x */
	AtVect x,	/* in: vector */
	double g,	/* in: multiplicand for y */
	AtVect y,	/* in: vector */
	AtVect z);	/* out: answer*/

/*
 * calc angular distance (RIKAKU)
 */
int atAngDistance(
	AtVect x,	/* input */
	AtVect y, 	/* input */
	double *r);	/* output: angular distance between x and y
				in radian */

/*
 * calc crossing points of cones (CROSS)
 */
int atCrossPts(
	AtVect x, 	/* input */
	double r1, 	/* input: angular distance from x (radian) */
	AtVect y, 	/* input */
	double r2, 	/* input: angular distance from y (radian) */
	AtVect z[2]);	/* output: two crossing points of the two cones */

/*
 * calc Vector Product (Gaiseki) of two vectors (SPOLE)
 */
int atVectProd(
	AtVect x, 	/* input */
	AtVect y, 	/* input */
	AtVect z);	/* output: vector (outer) product of x and y */

/*
 * calc Scalar Product (Naiseki, Dot-prooduct) of two vectors
 */
double atScalProd(	/* output: scaler (inner) product of x and y */
	AtVect x, 	/* input */
	AtVect y);	/* input */

/*
 * set Euler angles (EURST2)
 */
int atSetEuler(
	AtPolarVect *z, 	/* input: new Z-axis in old coordinate */
	AtPolarVect *y, 	/* input: vector in new Z-(+Y) plane*/
	AtEulerAng *ea);	/* output:Z-Y-Z Euler angle for the rotation */

/*
 * rotate polar vector with Euler angles (ROTAT2)
 */
int atRotPVect(
	AtEulerAng *ea, 	/* input */
	AtPolarVect *x, 	/* input */
	AtPolarVect *y);	/* output: result */

/* RECPOL */
/*
 * convert coordinate from Cartesian to polar (RECPOL)
 */
int atVectToPol(
	AtVect x, 		/* input */
	AtPolarVect *y);	/* output: result */

/*
 * convert coordinate from polar to Cartesian(POLREC)
 */
int atPolToVect(
	AtPolarVect *x, 	/* input */
	AtVect y);		/* output: result */

/*
 * convert R.A. and Dec from hh.mm.ss to radian
 */
int atConvPol(
	AtPolarVect60 *x, 	/* input */
	AtPolarVect *y); 	/* result */

/************************************************************************/
/*									*/
/*	Time related Functions						*/
/*									*/
/************************************************************************/

/*
 * convert AtTime to AtTimeD
 */
void atAtTimeToAtTimeD(
	AtTime *attime,		/* input: yr,mo,dy,hr,mn,sc,ms */
	AtTimeD *attimeD);	/* output: yr,mo,dy,hr,mn,sc,ss */

/*
 * convert AtTimeD to AtTime
 */
void atAtTimeDToAtTime(
	AtTimeD *attimeD,	/* input: yr,mo,dy,hr,mn,sc,ss */
	AtTime *attime);	/* input: yr,mo,dy,hr,mn,sc,ms */

/*
 * reformat AtTime, as 0<=ss<1, 0<=sc<=59, 0<=mn<=59, 0<=hr<=23, etc
 */
int atReformatAtTime(
	AtTime *time);	    /* input/output: yr,mo,dy,hr,mn,sc,ms */

int atReformatAtTimeD(
	AtTimeD *time);	    /* input/output: yr,mo,dy,hr,mn,sc,ss */

/*
 * format time in character (yy/mm/dd hh:mm:ss.ssssss) char[25]
 */
int atCTime(
	AtTime *time,		/* input: yr,mo,dy,hr,mn,sc,ms */
	char *ctime);		/* output: formated date (char [25]) */

int atCTimeD(
	AtTimeD *time,		/* input: yr,mo,dy,hr,mn,sc,ss */
	char *ctime);		/* output: formated date (char [25]) */

/*
 * format time in character (yyyy/mm/dd hh:mm:ss.ssssss) char[27]
 */
int atCTime2(
	AtTime *time,		/* input: yr,mo,dy,hr,mn,sc,ms */
	char *ctime);		/* output: formated date (char [27]) */

int atCTimeD2(
	AtTimeD *time,		/* input: yr,mo,dy,hr,mn,sc,ss */
	char *ctime);		/* output: formated date (char [27]) */

/* MJD */
/*
 * convert UT to Modified Julian Day (MJULIA)
 */
int atMJulian(
	AtTime *time,		/* input: yr,mo,dy,hr,mn,sc,ms */
	double *mjd);		/* output: modified julian day */

int atMJulianD(
	AtTimeD *time,		/* input: yr,mo,dy,hr,mn,sc,ss */
	double *mjd);		/* output: modified julian day */

/*
 * convert Modified Julian Day to UT (MJDATE)
 */
int atMJDate(
	double mjd,	    /* input: modified julian day */
	AtTime *time);	    /* output: yr,mo,dy,hr,mn,sc,ms */

int atMJDateD(
	double mjd,	    /* input: modified julian day */
	AtTimeD *time);	    /* output: yr,mo,dy,hr,mn,sc,ms */

/* Mission Time (elapsed seconds since specified MJD) */

/*
 * reset internal table
 */
void atMissionTimeResetTable(void);

/*
 * read leapsec table (or just return file name)
 */
char *atMissionTimeInit(	/* return: leapsec file name, NULL for ERROR */
	char *filename,	/* input: leapsec filename to read,
				  NULL for query,
				  "none" (case sensitive) for internal table */
	int flag);	/* input: +1: print leapsec table contents
				   0: do not print leapsec table
				  -1: print table, accept only FITS
				  -2: do not print table, accept only FITS */

/*
 * convert AtTimeD into mission time
 */
int atAtTimeDToMission(
	double mjd_base,	/* input: origin of mission time */
	AtTimeD *attime,	/* input: yr,mo,dy,hr,mn,sc,ss */
	double *mission);	/* output: mission time */

/*
 * convert mission time to AtTimeD
 */
int atMissionToAtTimeD(
	double mjd_base,	/* input: origin of mission time */
	double mission,	    	/* input: mission time */
	AtTimeD *attime);	/* output: yr,mo,dy,hr,mn,sc,ss */

/*
 * convert Modified Julian Day to mission time
 */
int atMJDToMission(
	double mjd_base,	/* input: origin of mission time */
	double mjd,		/* input: MJD */
	double *mission_time);	/* output: mission time */

/*
 * convert mission time to Modified Julain Day
 */
int atMissionToMJD(
	double mjd_base,	/* input: origin of mission time */
	double mission_time,	/* input: mission time */
	double *mjd);	    	/* output: MJD */

/************************************************************************/
/*									*/
/*	Functions related to ephemeis and orbit				*/
/*									*/
/************************************************************************/

/* definitions of coordinate used in the description of functions:

	sidereal: coordinate system fixed to the sky,
		refers to J2000 equatorial coordinate unless specified.

	geodetic: geocentric coordinate system  bound to the Earth.
		North pole is on Z-axis, (0,0) (equator at Greenwich longitude)
		defines X-axis, and y-axis is at east longitude = 90 deg.

	ground:	ground coordinate system at a tracking station on the Earth
	 	surface. Zenith defines the Z-axis, X-axis in the north,
		and Y in the west.

	azimuth and elevation: ground coordinate in polar representation, but
		the azimuth is measured from north as zero, then rotates
		in clockwise to the east, south, and to the west.

	geographical: Position of a place on the Earth in
		    "longitude", "latitude", and "altitude" in usual
		    geographical representation.

*/

/* SAISA */
/*
 * convert equatorial coordinate systems correcting for precession (SAISA)
 */
int atPrecession(
	double mjd0,	/* input: modified Julian Day */
	AtVect x0,	/* input: vector in equatorial coordiante at mjd0 */
	double mjd,	/* input: modified Julian Day */
	AtVect x);	/* output: vector in equatorial coordiante at mjd1 */

/*
 * Find Rotation Matrix for conversion of equatorial coordinate systems
 * correcting for precession (SAISA)
 */
int atPrecessRM(
	double mjd0,	/* input: modified Julian Day for the original coord*/
	double mjd,	/* input: modified Julian Day for the new coord*/
	AtRotMat rm);	/* output: rotation matrix to correct precession*/

/*
 * Find Euler Angles for conversion of equatorial coordinate systems
 * correcting for precession (SAISA)
 */
int atPrecessEuler(
	double mjd0,	/* input: modified Julian Day for the original coord*/
	double mjd,	/* input: modified Julian Day for the new coord*/
	AtEulerAng *ea);	/* output: Euler Ang to correct precession*/

/*
 * Equivalent to atPrecessRM(mjd, MJD_J2000, rm)
 */
void atPrecessRMJ2000(double mjd, AtRotMat rm);

/*
 * conversion of equatorial coordinate systems
 * correcting for precession (based on SAISA in ASTRO-C csub1)
 */

int atSaisa(
	double mjd0,		/* input: MJD for the original coord*/
	AtPolarVect *pv0,	/* input: original polar coordinate */
	double mjd,		/* input: MJD for the new coord*/
	AtPolarVect *pv);	/* output: Euler Ang to correct precession*/

/* ORBIT */
/*
 * set orbital elements (ELMST2): for plain text orbit file (OBSOLETE)
 */
int atSetElement(	/* return value: condition code */
	char *filename,	/* input: path name of the orbital element file */
	double mjd0,	/* input: modified Julian Day */
	int kchk);	/* input: if 0, calc 2nd derivative of mean anom. */

/*
 * set orbital elements: FITS version, link with "FITSIO" library  (OBSOLETE)
 */
int atSetElement2(	/* return value: condition code */
	char *filename,	/* input: FRF orbital element file name*/
	double mjd0,	/* input: modified Julian Day */
	int kchk);	/* input: if 0, calc 2nd derivative of mean anom. */

/*
 * calc path number (PATNUM) (OBSOLETE)
 */
int atPathNum(
	double mjd,	/* input: time in MJD */
	char path[11]);	/* output: path name */

/*
 * calc JAXA style path number "YYMMDDnnnn" at USC 34m station
 *	atSetElement() must be called before using this function (OBSOLETE)
 */
int atPathNumUSC(
	double mjd,	/* input: time in MJD */
	char path[11]);	/* output: path name, "YYMMDDnnnn" */

/*
 * calc satellite position (SATPNT) with the Earth center as origin (OBSOLETE)
 */
int atSatPos(
	double mjd,	/* input: time in MJD */
	AtVect x);	/* output: vector to the satellite */

/*
 * set orbital elements, with mjdz <= mjd0 < mjdn, see atElementTime3()
 */
int atSetElement3(	/* return value: condition code */
	AtElement3 *el,	/* input: contents of orbital elements to be read */
	char *filename,	/* input: FRF orbital element file name*/
	double mjd0,	/* input: modified Julian Day */
	int kchk);	/* input: if 0, calc 2nd derivative of mean anom. */

/*
 * get time of orbital elements, with mjdz <= mjd0 < mjdn, see atSetElement3()
 */
int atElementTime3(
	AtElement3 *el,	/* input: contents of orbital elements */
	double *mjdz,	/* output: current Epoch of orbital elements in MJD */
	double *mjdn);	/* output: next Epoch of orbital elements in MJD */

/*
 * calc path number (PATNUM)
 */
int atPathNum3(
	AtElement3 *el,	/* input: contents of orbital elements */
	double mjd,	/* input: time in MJD */
	char path[11]);	/* output: path name */

/*
 * calc JAXA style path number "YYMMDDnnnn" at USC 34m station
 */
int atPathNumUSC3(
	AtElement3 *el,	/* input: contents of orbital elements */
	double mjd,	/* input: time in MJD */
	char path[11]);	/* output: path name, "YYMMDDnnnn" */

/*
 * calc satellite position (SATPNT) with the Earth center as origin
 */
int atSatPos3(
	AtElement3 *el,	/* input: contents of orbital elements */
	double mjd,	/* input: time in MJD */
	AtVect x);	/* output: vector to the satellite */

/*
 * convert sidereal position to geodetic/geographical position
 * on earth (GEODCR)
 */
int atGeodcr(
	double mjd,	/* input: time in MJD */
	AtVect x,	/* input: vector in equatorial coordinate */
	double *heigh,	/* output: altitude [km] from the earth surface */
	double *longi,	/* output: longitude [rad] on the earth surface */
	double *latt);	/* output: latitude [rad] on the earth surface */

int atGeodetic(
	double mjd,	/* input: time in MJD */
	AtVect x,	/* input: vector in sidereal equatorial coordinate */
	AtVect y);	/* output: vector in geodetic coordinate at mjd */

int atInvGeodetic(
	double mjd,	/* input: time in MJD */
	AtVect x,	/* input: vector in geodetic coordinate at mjd */
	AtVect y);	/* output: vector in sidereal equatorial coordinate */

/*
 * convert polar geodetic coordinate latitude radial distance
 * to the geographic latitude and altitude from the earth surface
 * correcting for the ellipsoidal shape of the earth
 */
int atEllipsoid(
	AtPolarVect *xp, /* input: vector in celestial coordinate */
	double *latt,	/* output: latitude on the earth surface */
	double *heigh);	/* output: altitude from the earth surface */

/*
 * calc rotation matrix for converting equatorial position (J2000)
 *	to the geocentric position on earth
 */
int atSetGeoRM(
	double mjd,	/* input: time in MJD */
	AtRotMat rm);	/* output: rotation matrix */

/*
 * calc sidereal time (KOSEJI)
 */
int atSidereal(
	double mjd,		/* input: time in MJD */
	double *gsttod); 	/* output: Greenwich sidereal time (radian)
				 at mjd true of date*/

/*
 * solve Kepler equation (KEPLER)  g + e sin E = E
 */
int atKepler(
	double g,	/* input: mean anomaly */
	double eccent,	/* input: eccentricity */
	double *e);	/* output: eccentric anomaly */

/*
 * coordinate conversion from orbital plane to celestial (SATINA)
 */
int atOrbPlane(
	AtVect x,	/* input: vector in orbital plane */
	double ol,	/* input: Angle from Ascending Node */
	double ob, 	/* input: Right Ascention of the Ascending Node */
	double ai,	/* input: inclination of the orbit */
	AtVect y);	/* output: vector in celestial coordinate at mjd */

/*
 * convert geographic location of Tracking Station on Earth to Geodetic
 */
int atGeographic(
	AtPolarVect *y,	/* input: tracking station position;
				y.r:altitude from Earth surface
				y.lon, y.lat: longitude and lattitude */
	AtVect z);	/* output: geodetic coordinate of station */

/* same with atGeographic(), calling it internally */
int atGeographicToGeodetic(
	AtPolarVect *y,	/* input: tracking station position;
				y.r: altitude from Earth surface
				y.lon, y.lat: longitude and lattitude */
	AtVect z);	/* output: geodetic coordinate of station */

/* same with atGeographic(), calling it internally */
int atGeodesic(
	AtPolarVect *y,	/* input: tracking station position;
				y.r: altitude from Earth surface
				y.lon, y.lat: longitude and lattitude */
	AtVect z);	/* output: geodetic coordinate of station */

/*
 * convert geodetic vector into geographic position on earth surface
 */
int atGeodeticToGeographic(
	AtVect vect,	  /* input: geodetic coordinate of station */
	AtPolarVect *pv); /* output: geographic position on earth surface:
				pv->r: altitude from Earth surface
				pv->lon, pv->lat: longitude and lattitude */

/*
 * Convert from Geodetic to Ground Coordinate at the Tracking Station
 */
int atGroundCoord(
	AtVect station,	/* input: tracking station in geographic coord */
	AtRotMat rm);		/* output: rotation matrix to local coord*/

/*
 * Set Up geodetic vector to the tracking station and Rotation Matirix
 * to convert sidefloat coord to Ground (Earth-bound) coord
 *       Based on "ADAZEL"
 */

int atAzElSet(
	AtPolarVect *y,	/* input: geographic tracking station position :
				y.r:altitude from Earth surface
				y.lon, y.lat: longitude and lattitude */
	AtVect vStation,	/* output: geodetic coordinate of station */
	AtRotMat stationRM); /* output: rotation matrix which converts
				geographic coord to local coord at station*/
/*
 * Azimuth and Elevation for a vector in Ground (Earth-bound) coord
 */
int atAzEl(
	AtVect x,	    /* input: satellite pos. in geodetic coord. */
	AtVect stationV,    /* input: station pos. in geodetic coord. */
	AtRotMat stationRM, /* input: rotation matrix which converts geographic
				coord of tracking station to local coord */
	AtPolarVect *y);    /* output: local coordinate:
				y->r: distance, lon:azimuth, lat:elevation */

/*
 * examine earth occultation of specified direction (YEARTH)
 */
int atEarthOccult(
	AtVect satVect,	/* input: satellite pos. in equatorial coord. */
	AtVect xVect,	/* input: direction to be examined. */
	AtVect sunVect,	/* input: vector to the sun */
	int *flag, 	/* output: condition 0:sky, 1:dark, 2:bright earth */
	double *el);	/* output: elevation angle of z from the earth edge */

/* [OBSOLETE]
 * Set up Cut-off Rigidity Table (RIGSET), not required since atFunctions-2.7
 */
int atRigSet(
	char *filename);  /* input: path name of the cut-off rigidity file */

/* [OBSOLETE]
 * calc Cut-off Rigidity (RIGIDY), use atRigidityD() instead
 */
int atRigidity(
	AtPolarVect *x,	/* input: polar vector for geodetic position */
	float *rig);	/* output: cut-off rigidity (GeV/c) */

/* [OBSOLETE]
 * calc Cut-off Rigidity (RIGIDY)
 */
int atRigidityD(
	AtPolarVect *x,	/* input: polar vector for geodetic position */
	double *rig);	/* output: cut-off rigidity (GeV/c) */

/*
 * Set up Cut-off Rigidity Table in FITS format
 */
int atRigSet2(
	AtRigData2 **rdpp,	/* output: pointer to rigidity data pointer */
	char *filename);	/* input: FITS file name */

/*
 * Free allocated memory by atRigSet2()
 */
int atRigFree2(AtRigData2 *rdp);

/*
 * calc Cut-off Rigidity (RIGIDY)
 */
int atRigidity2(
	AtRigData2 *rdp,/* input: rigidity data pointer by atRigSet2() */
	AtPolarVect *x,	/* input: location [km] in geocentric coordinates */
	double *rig);	/* output: cut-off rigidity (GeV/c) */

/*
 * Set up Table of Geomagnetic Field (GASET)
 */
int atGeomagSet(
	double mjd,	/* input: time in MJD */
	int nmax);	/* input: maximum order of spherical harmonics */

/*
 * calc Geomagnetic Field (GA)
 */
int atGeomag(
	AtPolarVect *x,	/* input: polar vector for geodetic position */
	AtVect xs,	/* input: sidereal vector to the satellite */
	AtVect field);	/* output: sidereal direction of field line */

/*
 * calc position of planets, sun and moon (PLANETS, MOON in SOLARSYS)
 */
int atSun(
	double mjd,	/* input: time in MJD */
	AtVect pos);	/* output: vector to the sun in A.U.*/

int atPlanet(
	double mjd,	/* input: time in MJD */
	AtVect pos[9],	/* output: vector to the planets and sun (A.U.) */
	double size[9],	/* output: their visual size (radian)*/
	double mag[9]);	/* output: their visual magnitude*/

/*
 * calc position of the moon from the earth center in TOD coordinate.
 */
int atMoon(
	double mjd,	/* input: time in MJD */
	AtVect pos,	/* output: vector to the moon (km) */
	double *size,	/* output: their visual size (radian)*/
	double *phase,	/* output: phase (radian 0:new, pi:full)*/
	double *distan);/* output: distance to moon (km)*/

/*
 * calc if the point is in the "Brazil Anomaly" (SAA)
 */
int atSAA(
	int saa_type,	/* input: 0=nominal, 1=SIS, 2=STT, 3=HXD */
	AtPolarVect *x,	/* input: polar vector for geodetic position */
	int *flag);	/* output: =1 if in SAA, =0 if outside */

int atBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag);	/* output: =1 if in SAA, =0 if outside */

int atSISBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag);	/* output: =1 if in SAA, =0 if outside */

int atSTTBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag);	/* output: =1 if in SAA, =0 if outside */

int atHXDBrazil(
	double lon,	/* input: geodetic longitude in radian */
	double lat,	/* input: geodetic latitude in radian */
	int *flag);	/* output: =1 if in SAA, =0 if outside */

/************************************************************************/
/*									*/
/*	New for Astro-D 						*/
/*									*/
/************************************************************************/

/*
 * convert geocentric vector to Barycentric
 */
int atBarycentric(
	double mjd,	/* input: time in MJD */
	AtVect x,	/* input: vector in geocentric coordinate */
	AtVect y);	/* output: vector in barycentric coordinate at mjd */

/*
 * convert MJD to TDT (terrestial dynamical time)
 */
int atTDT(
	double mjd,	/* input: time in MJD */
	double *tdt);	/* output: time in TDT */

/************************************************************************/
/*									*/
/*	since version 1.8 						*/
/*									*/
/************************************************************************/

/*
 * examine (day|night) earth elevation of specified direction
 */
int atEarthElev(
	AtVect vSat,	/* input: satellite pos. in equatorial coord. */
	AtVect nvTgt,	/* input: direction to be examined. */
	AtVect vSun,	/* input: vector to the sun */
	int *occlt,	/* output: condition 0:sky, 1:dark, 2:bright earth */
	double elv[3]);	/* output: elevation [rad] from the earth/DYE/NTE */

/*
 * correct aberration
 */
void atAberration(
	double mjd,	/* input: modified Julian Day */
	AtVect x0,	/* input: vector in equatorial coordiante at mjd */
	AtVect x);	/* output: vector in equatorial coordiante at mjd */

void atInvAberration(
	double mjd,	/* input: modified Julian Day */
	AtVect x0,	/* input: vector in equatorial coordiante at mjd */
	AtVect x);	/* output: vector in equatorial coordiante at mjd */

/*
 * conversion between J2000, B1950, Galactic coordinates
 */
void atJ2000toB1950(double lo2000,double la2000,double *lo1950,double *la1950);
void atB1950toJ2000(double lo1950,double la1950,double *lo2000,double *la2000);
void atJ2000toGal(double ra, double dec, double *gl, double *gb);
void atGaltoJ2000(double gl, double gb, double *ra, double *dec);
void atGaltoB1950(double gl, double gb, double *ra, double *dec);
void atB1950toGal(double ra, double dec, double *gl, double *gb);

/*
 * conversion between Equatorial and Ecliptic coordinates
 */
void atJ2000ToEcliptic(double alpha, double delta, double *l, double *b);
void atEclipticToJ2000(double l, double b, double *alpha, double *delta);

/*
 * covert radian/degree to Right Ascension in hr,min,s.
 */
void atRadianToRA(
	double radian,		/* input: radian */
	AtRightAscension *ra);	/* output: Right Ascension */

void atDegToRA(
	double deg,		/* input: degree */
	AtRightAscension *ra);	/* output: Right Ascension */

/*
 * covert radian/degree to Declination in deg,min,s.
 */
void atRadianToDec(
	double radian,		/* input: radian */
	AtDeclination *dec);	/* output: Declination */

void atDegToDec(
	double deg,		/* input: degree */
	AtDeclination *dec);	/* output: Declination */

double atParseRAToRadian(	/* return: Right Ascension in radian */
	char *expression);	/* input: number (degree) or 00h00m00.0s */

double atParseDecToRadian(	/* return: Declination in radian */
	char *expression);	/* input: number (degree) or +00d00m00.0s */

/************************************************************************/
/*									*/
/*	since version 2.5 						*/
/*									*/
/************************************************************************/

/*
 * convert quaternion parameter set to Euler angles
 */
int atQuatToEuler(
	AtQuat q,		/* input: quaternion */
	AtEulerAng *ea);	/* output: z-y-z Euler Angle (radian) */

/*
 * convert Euler angles to quaternion parameter set
 */
int atEulerToQuat(
	AtEulerAng *ea,		/* input: z-y-z Euler Angle (radian) */
	AtQuat q);		/* output: quaternion */

/*
 * interpolation between two sets of q-parameters
 */
int
atInterpolateQuat(
	double t0, AtQuat q0,	/* input: q-parameter q0 at time t0 */
	double t1, AtQuat q1,	/* input: q-parameter q0 at time t0 */
	double t,  AtQuat q	/* input: time t
				   output: interpolated q-parameter q */
);

/*
 * interpolation between two sets of Euler angles
 */
int
atInterpolateEuler(
	double t0, AtEulerAng *ea0,	/* input: q-parameter q0 at time t0 */
	double t1, AtEulerAng *ea1,	/* input: q-parameter q0 at time t0 */
	double t,  AtEulerAng *ea	/* input: time t
					   output: interpolated Euler angles */
);


#ifdef __cplusplus
}
#endif

#endif
