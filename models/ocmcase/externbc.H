#ifndef _EXTERNBC_H_
#define _EXTERNBC_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <userfuncs.H>
#include <globalDefines.H>
#include <mflo_constants.H>
#include <species.H>

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
    s_ext[DENS_INDX] = s_int[DENS_INDX];
    s_ext[TEMP_INDX] = s_int[TEMP_INDX];
    s_ext[VELX_INDX] = s_int[VELX_INDX];
    s_ext[VELY_INDX] = s_int[VELY_INDX];
    s_ext[VELZ_INDX] = s_int[VELZ_INDX];

    s_ext[RHO_INDX]  = s_int[RHO_INDX];
    s_ext[RHOU_INDX] = s_int[RHOU_INDX];
    s_ext[RHOV_INDX] = s_int[RHOV_INDX];
    s_ext[RHOW_INDX] = s_int[RHOW_INDX];
    s_ext[RHOE_INDX] = s_int[RHOE_INDX];
    s_ext[VFRAC_INDX] = s_int[VFRAC_INDX];

    Real val;
    if(sgn > 0)
    {
        amrex::Real pres=101325.0; //Pa
        amrex::Real temp=600.0; //K
        amrex::Real Y_CH4=0.16;  //mass fractions
        amrex::Real Y_O2=0.08;
        amrex::Real Y_Ar=0.76;

        amrex::Real sum_YibyMi=0.0;
        amrex::Real R_gas;
        amrex::Real conc_CH4,conc_O2,conc_Ar;

        //P=rho (sum_i Yi/Mi) R T 

        sum_YibyMi += Y_CH4/mflo_species::molwts[CH4_ID];
        sum_YibyMi += Y_O2/mflo_species::molwts[O2_ID];
        sum_YibyMi += Y_Ar/mflo_species::molwts[Ar_ID];
        
        R_gas = sum_YibyMi*RU; //RU is universal gas constant

        amrex::Real rho=pres/temp/R_gas;

        conc_CH4 = rho*Y_CH4/mflo_species::molwts[CH4_ID];
        conc_O2  = rho*Y_O2/mflo_species::molwts[O2_ID];
        conc_Ar  = rho*Y_Ar/mflo_species::molwts[Ar_ID];


        //all species have 0 concentration including sites at this 
        //boundary
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
           s_ext[FLO_NVARS+sp]=0.0;
        }

        s_ext[FLO_NVARS+CH4_ID]=two*conc_CH4-s_int[FLO_NVARS+S1_ID];
        s_ext[FLO_NVARS+O2_ID]=two*conc_O2-s_int[FLO_NVARS+S1_ID];
        s_ext[FLO_NVARS+Ar_ID]=two*conc_Ar-s_int[FLO_NVARS+S1_ID];
    }
    else
    {
        //zero gradient
        for(int sp=0;sp<NUM_SPECIES;sp++)
        {
           s_ext[FLO_NVARS+sp]=s_int[FLO_NVARS+sp];
        }
    }
}
#endif
