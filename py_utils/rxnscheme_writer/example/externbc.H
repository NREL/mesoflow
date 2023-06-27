#ifndef _EXTERNBC_H_
#define _EXTERNBC_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <cflo_constants.H>
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
    //left/bottom/back boundary
    s_ext[PRES_INDX] = s_int[PRES_INDX];
    s_ext[TEMP_INDX] = cflo_user_funcs::fs_temp;
    s_ext[DENS_INDX] = cflo_thermo::get_r_from_tpc(cflo_user_funcs::fs_temp,
            s_ext[PRES_INDX],cflo_user_funcs::fs_spec);
    s_ext[VELX_INDX] = cflo_user_funcs::fs_vel;
    s_ext[VELY_INDX] = zero;
    s_ext[VELZ_INDX] = zero;
    s_ext[VFRAC_INDX] = one;


    s_ext[RHO_INDX]  = s_ext[DENS_INDX];
    s_ext[RHOU_INDX] = s_ext[DENS_INDX]*s_ext[VELX_INDX];
    s_ext[RHOV_INDX] = zero;
    s_ext[RHOW_INDX] = zero;
    
    s_ext[RHOE_INDX] = s_ext[DENS_INDX] *cflo_thermo::get_e_from_rpc(s_ext[DENS_INDX],
                                                  s_ext[PRES_INDX],
                                                  cflo_user_funcs::fs_spec);
    
    s_ext[RHOE_INDX]  += half*s_ext[DENS_INDX]*s_ext[VELX_INDX]*s_ext[VELX_INDX];

   

    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        s_ext[CFLO_NVARS+sp]=cflo_user_funcs::fs_spec[sp];
    }
}
#endif
