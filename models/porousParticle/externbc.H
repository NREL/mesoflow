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
    Real spec_ext[NUM_SPECIES];
 
    if(idir==1)
    {
        if(sgn > 0.0)
        {
            s_ext[DENS_INDX] = mflo_user_funcs::dens_left;
            s_ext[PRES_INDX] = mflo_user_funcs::pres_left;
            Real velmag=sqrt( s_int[VELX_INDX]*s_int[VELX_INDX]
                             + s_int[VELY_INDX]*s_int[VELY_INDX]
                             + s_int[VELZ_INDX]*s_int[VELZ_INDX] );
            s_ext[VELX_INDX] = velmag;
            s_ext[VELY_INDX] = zeroval;
            s_ext[VELZ_INDX] = zeroval;
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                spec_ext[sp]=mflo_user_funcs::spec_left[sp];
            }
            s_ext[TEMP_INDX] = mflo_thermo::get_t_from_rpc(s_ext[DENS_INDX],
                                                           s_ext[PRES_INDX],spec_ext);
        }
        else
        {
            s_ext[DENS_INDX] = s_int[DENS_INDX];
            s_ext[PRES_INDX] = mflo_user_funcs::pres_right;
            s_ext[VELX_INDX] = s_int[VELX_INDX];
            s_ext[VELY_INDX] = s_int[VELY_INDX];
            s_ext[VELZ_INDX] = s_int[VELZ_INDX];
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                spec_ext[sp]=s_int[FLO_NVARS+sp];
            }
            s_ext[TEMP_INDX] = mflo_thermo::get_t_from_rpc(s_ext[DENS_INDX],
                                                           s_ext[PRES_INDX],spec_ext);

        }

        s_ext[RHO_INDX]  = s_ext[DENS_INDX];
        s_ext[RHOU_INDX] = s_ext[DENS_INDX]*s_ext[VELX_INDX];
        s_ext[RHOV_INDX] = s_ext[DENS_INDX]*s_ext[VELY_INDX];
        s_ext[RHOW_INDX] = s_ext[DENS_INDX]*s_ext[VELZ_INDX];

        s_ext[RHOE_INDX]  = s_ext[DENS_INDX]*mflo_thermo::get_e_from_rpc(s_ext[DENS_INDX],
                                                            s_ext[PRES_INDX],spec_ext);
        s_ext[RHOE_INDX] += half*s_ext[DENS_INDX]*
        ( s_ext[VELX_INDX]*s_ext[VELX_INDX]
         + s_ext[VELY_INDX]*s_ext[VELY_INDX]
         + s_ext[VELZ_INDX]*s_ext[VELZ_INDX] );

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            s_ext[FLO_NVARS+sp]=spec_ext[sp];
        }

        s_ext[VFRAC_INDX] = s_int[VFRAC_INDX];
    }
}
#endif
