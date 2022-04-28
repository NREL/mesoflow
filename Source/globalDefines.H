#ifndef _GLOBAL_DEFINES_H_
#define _GLOBAL_DEFINES_H_
#include<species.H>

#define RHO_IND 0
#define RHOU_IND 1
#define RHOV_IND 2
#define RHOW_IND 3
#define RHOE_IND 4

#define DENS_IND 0
#define VELX_IND 1
#define VELY_IND 2
#define VELZ_IND 3
#define PRES_IND 4

#define NCVARS 5
//#define GAMA 1.4

#define RHO_INDX 0
#define RHOU_INDX 1
#define RHOV_INDX 2
#define RHOW_INDX 3
#define RHOE_INDX 4
#define DENS_INDX 5
#define VELX_INDX 6
#define VELY_INDX 7
#define VELZ_INDX 8
#define PRES_INDX 9
#define TEMP_INDX 10
#define VFRAC_INDX 11
#define FLO_NVARS 12
#define TOTAL_NVARS FLO_NVARS+NUM_SPECIES

#define VISC_INDX 0
#define THCOND_INDX 1

#define XDIR 0
#define YDIR 1
#define ZDIR 2

#define SMALLVAL 1e-18

#endif
