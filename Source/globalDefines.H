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
#define POT_INDX 12
#define EFLDX_INDX 13
#define EFLDY_INDX 14
#define EFLDZ_INDX 15
#define AX_INDX 16
#define AY_INDX 17
#define AZ_INDX 18
#define GAXX_INDX 19
#define GAXY_INDX 20
#define GAXZ_INDX 21
#define GAYX_INDX 22
#define GAYY_INDX 23
#define GAYZ_INDX 24
#define GAZX_INDX 25
#define GAZY_INDX 26
#define GAZZ_INDX 27
#define BX_INDX 28
#define BY_INDX 29
#define BZ_INDX 30
#define FLO_NVARS (BZ_INDX+1)
#define TOTAL_NVARS FLO_NVARS+NUM_SPECIES

#define VISC_INDX 0
#define THCOND_INDX 1

#define XDIR 0
#define YDIR 1
#define ZDIR 2

#define SMALLVAL 1e-18

//these functions are to get x,y,z directions
//given a sweep direction
//these are used for viscous flux calculation
//if sweep direction is 0, normal-x, t1-y, t2-z
//if sweep direction is 1, normal-y, t1-z, t2-x
//if sweep direction is 2, normal-z, t1-x, t2-y
#define GET_XDIR(a) (a*(a-1)/2     - 2*a*(a-2))
#define GET_YDIR(a) ((a-1)*(a-2)/2 + a*(a-1))
#define GET_ZDIR(a) ((a-1)*(a-2)   - a*(a-2))

#define PERBC 0
#define ROBINBC 1

#endif
