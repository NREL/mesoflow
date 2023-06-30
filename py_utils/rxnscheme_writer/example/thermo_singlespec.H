#ifndef _THERMO_H_
#define _THERMO_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <globalDefines.H>
#include <mflo_constants.H>
#include <species.H>

using namespace amrex;
namespace thermo
{
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_t_given_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real R_bg=RU/mflo_species::molwts[BG_ID];
            return(pres/R_bg/rho);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
        {
            Real gamma_bg  = 1.4;
            return(rhoe*(gamma_bg-one));
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
       Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
       {
            Real gamma_bg  = 1.4;
            Real gamma_mix=gamma_bg;
            return(gamma_mix);
       }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real gamma_bg  = 1.4;
            return( pres/rho/(gamma_bg-one) );
        }

}
#endif
