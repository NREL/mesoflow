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
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_r_from_c(Real spec[NUM_SPECIES],Real bgasconc=zeroval)
        {
            Real rho = zeroval;     
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                rho += spec[sp]*mflo_species::molwts[sp];   
            }
            return(rho);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_bgasconc_from_rc(Real rho,Real spec[NUM_SPECIES])
        {
            //this function is used only when there is background species
            return(0.0);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_t_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real sum_ci   = zeroval;
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci   += spec[sp];   
            }
            return(pres/(sum_ci*RU));
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
        {
            //not using coke and sites
            //these are only present in solid

            //include background gas also
            Real gamma_spec[NUM_GAS_SPECIES]={zeroval};
            Real cv_spec[NUM_GAS_SPECIES]={zeroval};
            Real T;

            gamma_spec[AIR_ID]  = 1.4;
            gamma_spec[A_ID]  = 1.2;
            gamma_spec[B_ID]  = 1.2;

            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_spec[sp]=RU/(gamma_spec[sp]-one); //J/K/mol
            }

            Real cv_times_spec=zeroval;
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_times_spec += cv_spec[sp]*spec[sp]; //J/K/m3
            }
            
            T=rhoe/cv_times_spec;
            Real p=zeroval;
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                p += spec[sp]*RU*T;   
            }

            return(p);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real gamma_spec[NUM_GAS_SPECIES]={zeroval};
            Real x_i[NUM_GAS_SPECIES]={zeroval};  //mole frac

            gamma_spec[AIR_ID] = 1.4;
            gamma_spec[A_ID]  = 1.2;
            gamma_spec[B_ID]  = 1.2;

            //find total number of moles
            Real sum_ci=zeroval;
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci += spec[sp];
            }

            //find molefraction
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                x_i[sp]=spec[sp]/sum_ci;
            }

            Real term2=zeroval;
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                term2 += x_i[sp]/(gamma_spec[sp]-one);
            }

            Real gama_mix=one+pow(term2,-one);

            return(gama_mix);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real gamma_spec[NUM_GAS_SPECIES]={zeroval};
            Real cv_spec[NUM_GAS_SPECIES]={zeroval};
            Real temp;

            gamma_spec[AIR_ID]  = 1.4;
            gamma_spec[A_ID]  = 1.2;
            gamma_spec[B_ID]  = 1.2;
            
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_spec[sp] = RU/(gamma_spec[sp]-one); //J/K/mol
            }

            temp=get_t_from_rpc(rho,pres,spec);

            Real sum_cvi_ci=zeroval;
            for(int sp=AIR_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_cvi_ci += cv_spec[sp]*spec[sp];
            }

            return( sum_cvi_ci*temp/rho );
        }

}
#endif
