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
        //modify if conjugate heat transfer 
        //included
        amrex::Real temp=300.0;
        return(temp);
    }
    
    AMREX_GPU_HOST_DEVICE AMREX_INLINE Real 
    get_solid_rhoe_from_t(Real temp)
    {
        //modify if conjugate 
        //heat transfer included
        amrex::Real Cv_solid=1.0;
        amrex::Real rho_solid=1.0;
        return(rho_solid*Cv_solid*temp);
    }
    extern AMREX_GPU_DEVICE_MANAGED Real R_univ;
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_r_from_c(Real spec[NUM_SPECIES],Real bgasconc=zeroval)
        {
            Real rho = zeroval;     
            for(int sp=0;sp<NUM_SPECIES;sp++)
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
            Real c_bg_gas = spec[N2_ID];
            return(c_bg_gas);
        }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_t_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real sum_ci=zeroval;
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                sum_ci += spec[sp];   
            }
            return(pres/(sum_ci*RU));
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
        {
            Real gamma_N2=1.4;
            Real gamma_He=1.6;

            Real cv_He=RU/(gamma_He-one); //J/K/mol
            Real cv_N2=RU/(gamma_N2-one); //J/K/mol

            Real T;

            T=rhoe/(cv_He*spec[HE_ID]+cv_N2*spec[N2_ID]);

            Real p=zeroval;
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                p += spec[sp]*RU*T;   
            }

            return(p);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
       Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
       {
            Real x_i[NUM_SPECIES];  //mole frac
            Real gamma_N2=1.4;
            Real gamma_He=1.6;

            Real sum_ci=zeroval;
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                sum_ci += spec[sp];
            }
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                x_i[sp]=spec[sp]/sum_ci;
            }

            Real term2 = pow(x_i[N2_ID]/(gamma_N2-one) + x_i[HE_ID]/(gamma_He-one),-one);
            Real gama_mix=one+term2;

            return(gama_mix);
       }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real x_i[NUM_SPECIES];  //mole frac
            Real temp;
            Real gamma_N2=1.4;
            Real gamma_He=1.6;
            
            Real cv_He=RU/(gamma_He-one)/mflo_species::molwts[HE_ID]; //J/K/kg
            Real cv_N2=RU/(gamma_N2-one)/mflo_species::molwts[N2_ID]; //J/K/kg

            Real sum_ci=zeroval;
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                sum_ci += spec[sp];
            }
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                x_i[sp]=spec[sp]/sum_ci;
            }

            temp=get_t_from_rpc(rho,pres,spec);
                
            return( temp*(cv_He*x_i[HE_ID]+cv_N2*x_i[N2_ID]) );
        }

}
#endif
