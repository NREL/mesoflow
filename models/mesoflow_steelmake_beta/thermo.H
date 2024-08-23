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
get_solid_rhoe_from_t(Real temp)
{
    amrex::Real Cv_solid=820.0;
    amrex::Real rho_solid=5240.0;
    return(rho_solid*Cv_solid*temp);
}






// Calculate the density from species concentrations
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_r_from_c(Real spec[NUM_SPECIES],Real bgasconc=zeroval)
        {
            Real rho = zeroval;     
            for(int sp=0;sp<NUM_SPECIES;sp++)
            {
                //density = species concentrations * molecular weights
                rho += spec[sp]*mflo_species::molwts[sp];   
            }
            return(rho);
        }

// // Get the concentration of background gas from the density and species concentrations
//    AMREX_GPU_HOST_DEVICE AMREX_INLINE
//        Real get_bgasconc_from_rc(Real rho,Real spec[NUM_SPECIES])
//        {
//            //this function is used only when there is background species
//
//            return(0.0);
//        }


//
AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_bgasconc_from_rc(Real rho,Real spec[NUM_SPECIES])
    {
        //this function is used only when there is background species
        //that is not being tracked
        Real c_bg_gas = spec[N2_ID];
        return(c_bg_gas);
    }



// Calculate the temperature from density, pressure, and species concentrations
AMREX_GPU_HOST_DEVICE AMREX_INLINE
  Real get_t_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
  {
    Real sum_ci   = zeroval;
    for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
    {
        sum_ci   += spec[sp];
    }
    return(pres/(sum_ci*RU));
  }
 



// Calculate the specific heat ratio (gamma) from density, pressure, and species concentrations
AMREX_GPU_HOST_DEVICE AMREX_INLINE
 Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
 {
     Real gamma_spec[NUM_GAS_SPECIES]={zeroval};
     Real x_i[NUM_GAS_SPECIES]={zeroval};  //mole frac

     gamma_spec[H2_ID]  = 1.4;
     gamma_spec[H2O_ID]  = 1.3;
     gamma_spec[CO_ID]  = 1.4 ;
     gamma_spec[N2_ID]  = 1.4;
     gamma_spec[CO2_ID]  = 1.3;
    

     //find total number of moles
     Real sum_ci=zeroval;
     for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
     {
         sum_ci += spec[sp];
     }

     //find molefraction
     for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
     {
         x_i[sp]=spec[sp]/sum_ci;
     }

     Real term2=zeroval;
     for(int sp=0;sp<NUM_GAS_SPECIES;sp++)
     {
         term2 += x_i[sp]/(gamma_spec[sp]-one);
     }

     Real gamma_mix=one+pow(term2,-one);

     return(gamma_mix);
     
 }

// // Calculate the pressure from internal energy, density, and species concentrations
//    AMREX_GPU_HOST_DEVICE AMREX_INLINE
//        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
//        {
//            return(rhoe*(gamma_mix-one));
//        }



// // Calculate the internal energy from density, pressure, and species concentrations
//    AMREX_GPU_HOST_DEVICE AMREX_INLINE
//        Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
//        {
//            return( pres/(gamma_mix-one)/rho );
//        }




AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
    {
        Real gamma_N2=1.4;
        Real gamma_H2=1.4;
        Real gamma_H2O=1.3;
        Real gamma_CO=1.4;
        Real gamma_CO2=1.3;
        
        


        Real cv_H2=RU/(gamma_H2-one); //J/K/mol
        Real cv_N2=RU/(gamma_N2-one); //J/K/mol
        Real cv_H2O=RU/(gamma_H2O-one); //J/K/mol
        Real cv_CO=RU/(gamma_CO-one); //J/K/mol
        Real cv_CO2=RU/(gamma_CO2-one); //J/K/mol

        Real T;

        T=rhoe/(cv_H2*spec[H2_ID]+cv_CO*spec[CO_ID]+cv_H2O*spec[H2O_ID]+cv_CO2*spec[CO2_ID])+cv_N2*spec[N2_ID];

        Real p=zeroval;
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
            p += spec[sp]*RU*T;
        }

        return(p);
    }



AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
    {
        Real x_i[NUM_SPECIES];  //mole frac
        Real temp;
        Real gamma_N2=1.4;
        Real gamma_H2=1.4;
        Real gamma_H2O=1.3;
        Real gamma_CO=1.4;
        Real gamma_CO2=1.3;
       
        
        Real cv_N2=RU/(gamma_N2-one)/mflo_species::molwts[N2_ID]; //J/K/kg
        Real cv_H2=RU/(gamma_H2-one)/mflo_species::molwts[H2_ID]; //J/K/kg
        Real cv_CO=RU/(gamma_CO-one)/mflo_species::molwts[CO_ID]; //J/K/kg
        Real cv_CO2=RU/(gamma_CO2-one)/mflo_species::molwts[CO2_ID];
        Real cv_H2O=RU/(gamma_H2O-one)/mflo_species::molwts[H2O_ID];
        

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
            
        return( temp*(cv_H2*x_i[H2_ID]+ cv_CO*x_i[CO_ID]  + cv_H2O*x_i[H2O_ID]+cv_CO2*x_i[CO2_ID]+cv_N2*x_i[N2_ID] ) );//
    }






}
#endif
