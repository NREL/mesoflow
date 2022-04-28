#ifndef _EXTERNBC_H_
#define _EXTERNBC_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>

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
    s_ext[DENS_INDX] = mflo_user_funcs::rho_l;
    s_ext[TEMP_INDX] = mflo_user_funcs::T_l;
    s_ext[VELX_INDX] = zeroval;
    s_ext[VELY_INDX] = zeroval;
    s_ext[VELZ_INDX] = zeroval;
    s_ext[VFRAC_INDX] = one;

    s_ext[RHO_INDX]  = s_ext[DENS_INDX];
    s_ext[RHOU_INDX] = zeroval;
    s_ext[RHOV_INDX] = zeroval;
    s_ext[RHOW_INDX] = zeroval;
    s_ext[RHOE_INDX] = s_ext[DENS_INDX]*mflo_thermo::get_e_from_rpc(mflo_user_funcs::rho_l,
                                                               mflo_user_funcs::p_l,
                                                               mflo_user_funcs::spec_l);

    for(int sp=0;sp<NUM_SPECIES;sp++)
    {
        s_ext[FLO_NVARS+sp]=mflo_user_funcs::spec_l[sp];
    }
}
#endif
