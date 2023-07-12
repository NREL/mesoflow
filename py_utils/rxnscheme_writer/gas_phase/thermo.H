
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
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                rho += spec[sp]*mflo_species::molwts[sp];   
            }
            rho += bgasconc*BG_GAS_MWT;

            return(rho);
        }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_bgasconc_from_rc(Real rho,Real spec[NUM_SPECIES])
        {
            Real sum_ciMi=0.0;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*mflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;
            return(c_bg_gas);
        }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_r_from_tpc(Real temp,Real pres,Real spec[NUM_SPECIES])
        {
            Real sum_ci   = zeroval;
            Real sum_ciMi = zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci   += spec[sp];   
                sum_ciMi += spec[sp]*mflo_species::molwts[sp];   
            }

            Real c_total = pres/RU/temp;
            Real c_bg_gas = c_total-sum_ci;

            return(sum_ciMi+c_bg_gas*BG_GAS_MWT);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_t_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real sum_ci   = zeroval;
            Real sum_ciMi = zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci   += spec[sp];   
                sum_ciMi += spec[sp]*mflo_species::molwts[sp];   
            }

            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;
            sum_ci += c_bg_gas;

            //return(pres*BG_GAS_MWT/rho/RU);
            return(pres/(sum_ci*RU));
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_p_from_rec(Real rhoe,Real rho,Real spec[NUM_SPECIES])
        {
            //not using coke and sites
            //these are only present in solid

            //include background gas also
            Real cv_spec[NUM_GAS_SPECIES+1]={zeroval};
            Real T;

            //find background gas concentration
            Real sum_ciMi = zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*mflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;

            for(int sp=FIRST_SPEC;sp<(NUM_GAS_SPECIES+1);sp++)
            {
                cv_spec[sp]=RU/(mflo_species::gamma_spec[sp]-one); //J/K/mol
            }

            Real cv_times_spec=zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_times_spec += cv_spec[sp]*spec[sp]; //J/K/m3
            }
            cv_times_spec += cv_spec[NUM_GAS_SPECIES]*c_bg_gas;


            T=rhoe/cv_times_spec;
            Real p=zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                p += spec[sp]*RU*T;   
            }
            p += c_bg_gas*RU*T;

            //return(rhoe*(GAMMA_BG_GAS-one));
            return(p);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_gama_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real x_i[NUM_GAS_SPECIES+1]={zeroval};  //mole frac

            //find background gas concentration
            Real sum_ciMi = zeroval;
            for(int sp=Ar_ID;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*mflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;


            Real sum_ci=zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ci += spec[sp];
            }
            sum_ci += c_bg_gas;

            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                x_i[sp]=spec[sp]/sum_ci;
            }
            x_i[NUM_GAS_SPECIES] = c_bg_gas/sum_ci;

            Real term2=zeroval;
            for(int sp=FIRST_SPEC;sp<(NUM_GAS_SPECIES+1);sp++)
            {
                term2 += x_i[sp]/(mflo_species::gamma_spec[sp]-one);
            }

            Real gama_mix=one+pow(term2,-one);

            //return(GAMMA_BG_GAS);
            return(gama_mix);
        }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real get_e_from_rpc(Real rho,Real pres,Real spec[NUM_SPECIES])
        {
            Real cv_spec[NUM_GAS_SPECIES+1]={zeroval};
            Real temp;

            //find background gas concentration
            Real sum_ciMi = zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_ciMi += spec[sp]*mflo_species::molwts[sp];   
            }
            Real c_bg_gas = (rho-sum_ciMi)/BG_GAS_MWT;

            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                cv_spec[sp] = RU/(mflo_species::gamma_spec[sp]-one); //J/K/mol
            }
            cv_spec[NUM_GAS_SPECIES] = RU/(mflo_species::gamma_spec[NUM_GAS_SPECIES]-one); 

            temp=get_t_from_rpc(rho,pres,spec);

            Real sum_cvi_ci=zeroval;
            for(int sp=FIRST_SPEC;sp<NUM_GAS_SPECIES;sp++)
            {
                sum_cvi_ci += cv_spec[sp]*spec[sp];
            }
            sum_cvi_ci += cv_spec[NUM_GAS_SPECIES]*c_bg_gas;

            //return(pres/rho/(GAMMA_BG_GAS-one));
            return( sum_cvi_ci*temp/rho );
        }

}
#endif

        