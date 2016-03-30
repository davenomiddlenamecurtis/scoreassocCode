#ifndef SAFILTERFUNCSHPP
#define SAFILTERFUNCSHPP

#include "scoreassoc.hpp"

#ifndef USEFILTERS
#define USEFILTERS
#endif

int initExclusions(FILE *fp,char *extras[]=0);
int applyExclusions(par_info *pi);
int stateExclusions(FILE *fp);

extern subject **global_sub;
#endif