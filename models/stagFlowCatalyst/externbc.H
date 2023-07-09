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
    if( (x[YDIR]*x[YDIR]+x[ZDIR]*x[ZDIR]) 
       < mflo_user_funcs::jetrad*mflo_user_funcs::jetrad)
    {
        Real spec[NUM_SPECIES]={zeroval};
        s_ext[DENS_INDX]=mflo_user_funcs::rho0;
        s_ext[PRES_INDX]=s_int[PRES_INDX];
        s_ext[VELX_INDX]=mflo_user_funcs::fs_vel;
        s_ext[VELY_INDX]=zeroval;
        s_ext[VELZ_INDX]=zeroval;
        
        spec[CH4_ID] = mflo_user_funcs::CH4conc;
        spec[O2_ID] = (s_ext[DENS_INDX]
                       -mflo_user_funcs::CH4conc*mflo_species::molwts[CH4_ID])*0.23/mflo_species::molwts[O2_ID];
        spec[N2_ID] = (s_ext[DENS_INDX]
                       -mflo_user_funcs::CH4conc*mflo_species::molwts[CH4_ID])*0.77/mflo_species::molwts[N2_ID];
        spec[C2H4_ID] = zeroval;
        spec[H2O_ID]  = zeroval;
        spec[S_ID]   = zeroval;

        s_ext[TEMP_INDX] = mflo_thermo::get_t_from_rpc(s_ext[DENS_INDX],
                                                       s_ext[PRES_INDX],spec);

        s_ext[RHO_INDX]  = s_ext[DENS_INDX];
        s_ext[RHOU_INDX] = s_ext[DENS_INDX]*s_ext[VELX_INDX];
        s_ext[RHOV_INDX] = s_ext[DENS_INDX]*s_ext[VELY_INDX];
        s_ext[RHOW_INDX] = s_ext[DENS_INDX]*s_ext[VELZ_INDX];
        s_ext[RHOE_INDX]  = s_ext[DENS_INDX]*
        mflo_thermo::get_e_from_rpc(s_ext[DENS_INDX],s_ext[PRES_INDX],spec);

        s_ext[RHOE_INDX] += half*s_ext[DENS_INDX]*
        (s_ext[VELX_INDX]*s_ext[VELX_INDX] +
         s_ext[VELY_INDX]*s_ext[VELY_INDX] +
         s_ext[VELZ_INDX]*s_ext[VELZ_INDX]);

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            s_ext[FLO_NVARS+sp]=spec[sp];
        }
        s_ext[VFRAC_INDX] = s_int[VFRAC_INDX];
    }
    else
    {
        //wall bc
        s_ext[RHO_INDX]  = s_int[RHO_INDX];
        s_ext[RHOU_INDX] = -s_int[VELX_INDX];
        s_ext[RHOV_INDX] = -s_int[VELY_INDX];
        s_ext[RHOW_INDX] = -s_int[VELZ_INDX];
        s_ext[RHOE_INDX]  = s_int[RHOE_INDX];

        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            s_ext[FLO_NVARS+sp]=s_int[FLO_NVARS+sp];
        }
        s_ext[VFRAC_INDX] = s_int[VFRAC_INDX];

        s_ext[DENS_INDX]=s_int[DENS_INDX];
        s_ext[PRES_INDX]=s_int[PRES_INDX];
        s_ext[VELX_INDX]=-s_int[VELX_INDX];
        s_ext[VELY_INDX]=-s_int[VELY_INDX];
        s_ext[VELZ_INDX]=-s_int[VELZ_INDX];
    }
}
#endif
