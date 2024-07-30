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
namespace mflo_thermo
{
    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real 
    get_solid_t_from_rhoe(Real rhoe)
    {
        amrex::Real Cv_solid=1000.0;
        amrex::Real rho_solid=2000.0;
        amrex::Real temp=rhoe/(rho_solid*Cv_solid);
        return(temp);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_r_from_c(Real spec[NUM_SPECIES],Real bgasconc=zeroval)
        {
            Real rho = zeroval;     
            for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
            {
                rho += spec[sp]*mflo_species::molwts[sp];   
            }
            return(rho);
        }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_bgasconc_from_rc(Real rho,Real spec[NUM_SPECIES])
        {
            //this function is used only when there is background species
            //that is not being tracked
            Real c_bg_gas = spec[S1_ID];
            return(c_bg_gas);
        }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_t_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real Rgas=one;
            return(pres/rho/Rgas);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
        {
            Real gamma=1.4;
            return(rhoe*(gamma-one));
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
       Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
       {
            Real gamma=1.4;
            return(gamma);
       }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real gamma=1.4;
            return( pres/(gamma-one)/rho );
        }

}
#endif
