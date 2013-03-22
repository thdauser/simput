#include "atFunctions.h"
#include "atError.h"

int                   /* Giving Euler-angles to the new-coordinate (EURST2) */
atSetEuler(           /*                             ver 1.0  92/07/01  ay  */
        AtPolarVect *z,         /* input: new Z-axis in old coordinate */
        AtPolarVect *y,         /* input: vector in new Z-(+Y) plane*/
        AtEulerAng *ea)         /* output:Z-Y-Z Euler angle for the rotation */
{
  AtPolarVect v;
  int  rc;
  ea->phi = z->lon;  ea->theta = PI/2.0 - z->lat;  ea->psi = 0.0;
  if ( (rc=atRotPVect( ea, y, &v )) != 0 ) { return (rc); }
  ea->psi = v.lon - PI/2.0;
  return (NORMAL_END);
}
