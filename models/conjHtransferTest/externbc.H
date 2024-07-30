#ifndef _EXTERNBC_H_
#define _EXTERNBC_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>
#include <species.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void externalbc(const amrex::Real x[AMREX_SPACEDIM],
        const amrex::Real s_int[],
        amrex::Real s_ext[],
        const int idir,
        const int sgn,
        const amrex::Real time,
        amrex::GeometryData const& geomdata)
{
    //default adiabatic wall
    s_ext[PRES_INDX] = s_int[PRES_INDX];
    s_ext[DENS_INDX] = s_int[DENS_INDX];
    s_ext[TEMP_INDX] = s_int[TEMP_INDX];
    s_ext[VELX_INDX] = -s_int[VELX_INDX];
    s_ext[VELY_INDX] = -s_int[VELY_INDX];
    s_ext[VELZ_INDX] = -s_int[VELZ_INDX];

    s_ext[RHO_INDX]  = s_int[RHO_INDX];
    s_ext[RHOU_INDX] = -s_int[RHOU_INDX];
    s_ext[RHOV_INDX] = -s_int[RHOV_INDX];
    s_ext[RHOW_INDX] = -s_int[RHOW_INDX];
    s_ext[RHOE_INDX] = s_int[RHOE_INDX];
    s_ext[VFRAC_INDX] = s_int[VFRAC_INDX];
    
    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        s_ext[FLO_NVARS+sp]=s_int[FLO_NVARS+sp];
    }

    if(sgn > 0)
    {
        s_ext[TEMP_INDX]=1100.0;
    }
    else
    {
        s_ext[TEMP_INDX]=300.0;
    }
    
}
#endif
