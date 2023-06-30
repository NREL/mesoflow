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
            Real sum_ci=zeroval;
            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci += spec[sp];   
            }
            return(pres/(sum_ci*RU));
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
        {
            //not using coke and sites
            //these are only present in solid

            Real gamma_spec[NUM_GAS_SPECIES]={zeroval};
            Real cv_spec[NUM_GAS_SPECIES]={zeroval};
            Real T;

            gamma_spec[BG_ID]  = 1.4;
            gamma_spec[PV_ID]  = 1.2;
            gamma_spec[HC_ID]  = 1.2;
            gamma_spec[FPN_ID] = 1.2;

            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_spec[sp]=RU/(gamma_spec[sp]-one); //J/K/mol
            }

            Real cv_times_spec=zeroval;
            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_times_spec += cv_spec[sp]*spec[sp]; //J/K/m3
            }

            T=rhoe/cv_times_spec;

            Real p=zeroval;
            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
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

            gamma_spec[BG_ID]  = 1.4;
            gamma_spec[PV_ID]  = 1.2;
            gamma_spec[HC_ID]  = 1.2;
            gamma_spec[FPN_ID] = 1.2;

            Real sum_ci=zeroval;
            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci += spec[sp];
            }
            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                x_i[sp]=spec[sp]/sum_ci;
            }
            
            Real term2=zeroval;
            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
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
            Real x_i[NUM_GAS_SPECIES]={zeroval};  //mole frac
            Real cv_spec[NUM_GAS_SPECIES]={zeroval};
            Real temp;

            gamma_spec[BG_ID]  = 1.4;
            gamma_spec[PV_ID]  = 1.2;
            gamma_spec[HC_ID]  = 1.2;
            gamma_spec[FPN_ID] = 1.2;
           
            for(int sp=BG_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_spec[sp] = RU/(gamma_spec[sp]-one)/mflo_species::molwts[sp]; //J/K/kg
            } 

            Real sum_ci=zeroval;
            for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci += spec[sp];
            }
            for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
            {
                x_i[sp]=spec[sp]/sum_ci;
            }

            temp=get_t_given_rpc(rho,pres,spec);

            Real cvbar=zeroval;
            for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
            {
                cvbar += cv_spec[sp]*x_i[sp];
            }

            return( cvbar*temp );
        }

}
#endif
