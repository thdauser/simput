#ifndef _POSSTRING_H
#define _POSSTRING_H

#define DEGPOS  0
#define LATPOS  1
#define LONPOS  2
#define HOURPOS 3


double parse_posstring(const char *str, const int type);
char *posstring(double position, const int type);

#endif

